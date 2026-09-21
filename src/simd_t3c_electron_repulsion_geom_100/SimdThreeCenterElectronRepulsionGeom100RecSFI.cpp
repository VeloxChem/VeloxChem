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


#include "SimdThreeCenterElectronRepulsionGeom100RecSFI.hpp"

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
compute_geom_100_sfi_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_geom_100_sfi_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 22614, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 273 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 22614, 21644, 840, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 9, 6, 10,
                                                             ncols, fj, i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 21, 3, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 24, 3, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 27, 3, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 30, 3, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 33, 3, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 36, 3, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 39, 3, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 42, 3, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 45, 3, 6, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 48, 3, 6, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 51, 3, 6, 10, 11,
                                                                       21, 24, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 57, 3, 6, 11, 12,
                                                                       24, 27, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 63, 3, 6, 12, 13,
                                                                       27, 30, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 69, 3, 6, 13, 14,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 75, 3, 6, 14, 15,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 81, 3, 6, 15, 16,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 87, 3, 6, 16, 17,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 93, 3, 6, 17, 18,
                                                                       42, 45, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 99, 3, 6, 18, 19,
                                                                       45, 48, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 105, 3, 6, 21, 24,
                                                                       51, 57, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 115, 3, 6, 24, 27,
                                                                       57, 63, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 125, 3, 6, 27, 30,
                                                                       63, 69, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 135, 3, 6, 30, 33,
                                                                       69, 75, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 145, 3, 6, 33, 36,
                                                                       75, 81, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 155, 3, 6, 36, 39,
                                                                       81, 87, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 165, 3, 6, 39, 42,
                                                                       87, 93, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 175, 3, 6, 42, 45,
                                                                       93, 99, ncols, gamma, p,
                                                                       q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 185, 0, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 188, 0, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 191, 0, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 194, 0, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 197, 0, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 200, 0, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 203, 0, 6, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 206, 0, 6, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 209, 0, 6, 10, 11,
                                                                       21, 24, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 218, 0, 6, 11, 12,
                                                                       24, 27, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 227, 0, 6, 12, 13,
                                                                       27, 30, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 236, 0, 6, 13, 14,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 245, 0, 6, 14, 15,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 254, 0, 6, 15, 16,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 263, 0, 6, 16, 17,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 272, 0, 6, 17, 18,
                                                                       42, 45, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 281, 0, 6, 18, 19,
                                                                       45, 48, ncols, gamma, p,
                                                                       q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 290, 0, 3, 6, 21,
                                                                       24, 51, 57, 209, 218,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 308, 0, 3, 6, 24,
                                                                       27, 57, 63, 218, 227,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 326, 0, 3, 6, 27,
                                                                       30, 63, 69, 227, 236,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 344, 0, 3, 6, 30,
                                                                       33, 69, 75, 236, 245,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 362, 0, 3, 6, 33,
                                                                       36, 75, 81, 245, 254,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 380, 0, 3, 6, 36,
                                                                       39, 81, 87, 254, 263,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 398, 0, 3, 6, 39,
                                                                       42, 87, 93, 263, 272,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 416, 0, 3, 6, 42,
                                                                       45, 93, 99, 272, 281,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 434, 0, 3, 6, 51,
                                                                       57, 105, 115, 290, 308,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 464, 0, 3, 6, 57,
                                                                       63, 115, 125, 308, 326,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 494, 0, 3, 6, 63,
                                                                       69, 125, 135, 326, 344,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 524, 0, 3, 6, 69,
                                                                       75, 135, 145, 344, 362,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 554, 0, 3, 6, 75,
                                                                       81, 145, 155, 362, 380,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 584, 0, 3, 6, 81,
                                                                       87, 155, 165, 380, 398,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 614, 0, 3, 6, 87,
                                                                       93, 165, 175, 398, 416,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 644, 6, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 647, 6, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 650, 6, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 653, 6, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 656, 6, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 659, 6, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 662, 6, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 665, 6, 19, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 668, 6, 20, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 671, 6, 12, 27,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 680, 6, 13, 30,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 689, 6, 14, 33,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 698, 6, 15, 36,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 707, 6, 16, 39,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 716, 6, 17, 42,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 725, 6, 18, 45,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 734, 6, 19, 48,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 743, 6, 27, 63,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 761, 6, 30, 69,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 779, 6, 33, 75,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 797, 6, 36, 81,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 815, 6, 39, 87,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 833, 6, 42, 93,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 851, 6, 45, 99,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 869, 6, 63, 125,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 899, 6, 69, 135,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 929, 6, 75, 145,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 959, 6, 81, 155,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 989, 6, 87, 165,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1019, 6, 93, 175,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1049, 6, 12, 185,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1058, 6, 13, 188,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1067, 6, 14, 191,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1076, 6, 15, 194,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1085, 6, 16, 197,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1094, 6, 17, 200,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1103, 6, 18, 203,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1112, 6, 19, 206,
                                                                       ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1121, 6, 27, 185,
                                                                       227, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1148, 6, 30, 188,
                                                                       236, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1175, 6, 33, 191,
                                                                       245, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1202, 6, 36, 194,
                                                                       254, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1229, 6, 39, 197,
                                                                       263, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1256, 6, 42, 200,
                                                                       272, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1283, 6, 45, 203,
                                                                       281, ncols, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1310, 3, 6, 63,
                                                                       1121, 227, 1148, 326,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1364, 3, 6, 69,
                                                                       1148, 236, 1175, 344,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1418, 3, 6, 75,
                                                                       1175, 245, 1202, 362,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1472, 3, 6, 81,
                                                                       1202, 254, 1229, 380,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1526, 3, 6, 87,
                                                                       1229, 263, 1256, 398,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1580, 3, 6, 93,
                                                                       1256, 272, 1283, 416,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 1634, 3, 6, 125,
                                                                       1310, 326, 1364, 494,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 1724, 3, 6, 135,
                                                                       1364, 344, 1418, 524,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 1814, 3, 6, 145,
                                                                       1418, 362, 1472, 554,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 1904, 3, 6, 155,
                                                                       1472, 380, 1526, 584,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 1994, 3, 6, 165,
                                                                       1526, 398, 1580, 614,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2084, 6, 10, 11,
                                                                       644, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2090, 6, 11, 12,
                                                                       647, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2096, 6, 12, 13,
                                                                       650, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2102, 6, 13, 14,
                                                                       653, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2108, 6, 14, 15,
                                                                       656, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2114, 6, 15, 16,
                                                                       659, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2120, 6, 16, 17,
                                                                       662, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2126, 6, 17, 18,
                                                                       665, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2132, 6, 18, 19,
                                                                       668, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2138, 3, 6, 2084,
                                                                       644, 2090, 671, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2156, 3, 6, 2090,
                                                                       647, 2096, 680, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2174, 3, 6, 2096,
                                                                       650, 2102, 689, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2192, 3, 6, 2102,
                                                                       653, 2108, 698, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2210, 3, 6, 2108,
                                                                       656, 2114, 707, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2228, 3, 6, 2114,
                                                                       659, 2120, 716, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2246, 3, 6, 2120,
                                                                       662, 2126, 725, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2264, 3, 6, 2126,
                                                                       665, 2132, 734, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2282, 3, 6, 2138,
                                                                       671, 2156, 51, 57, 743,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2318, 3, 6, 2156,
                                                                       680, 2174, 57, 63, 761,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2354, 3, 6, 2174,
                                                                       689, 2192, 63, 69, 779,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2390, 3, 6, 2192,
                                                                       698, 2210, 69, 75, 797,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2426, 3, 6, 2210,
                                                                       707, 2228, 75, 81, 815,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2462, 3, 6, 2228,
                                                                       716, 2246, 81, 87, 833,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2498, 3, 6, 2246,
                                                                       725, 2264, 87, 93, 851,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2534, 3, 6, 2282,
                                                                       743, 2318, 105, 115, 869,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2594, 3, 6, 2318,
                                                                       761, 2354, 115, 125, 899,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2654, 3, 6, 2354,
                                                                       779, 2390, 125, 135, 929,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2714, 3, 6, 2390,
                                                                       797, 2426, 135, 145, 959,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2774, 3, 6, 2426,
                                                                       815, 2462, 145, 155, 989,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2834, 3, 6, 2462,
                                                                       833, 2498, 155, 165, 1019,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2894, 0, 6, 2084,
                                                                       644, 2090, 1049, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2912, 0, 6, 2090,
                                                                       647, 2096, 1058, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2930, 0, 6, 2096,
                                                                       650, 2102, 1067, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2948, 0, 6, 2102,
                                                                       653, 2108, 1076, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2966, 0, 6, 2108,
                                                                       656, 2114, 1085, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2984, 0, 6, 2114,
                                                                       659, 2120, 1094, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3002, 0, 6, 2120,
                                                                       662, 2126, 1103, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3020, 0, 6, 2126,
                                                                       665, 2132, 1112, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3038, 0, 3, 6,
                                                                       2138, 671, 2156, 2894,
                                                                       1049, 2912, 209, 218,
                                                                       1121, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3092, 0, 3, 6,
                                                                       2156, 680, 2174, 2912,
                                                                       1058, 2930, 218, 227,
                                                                       1148, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3146, 0, 3, 6,
                                                                       2174, 689, 2192, 2930,
                                                                       1067, 2948, 227, 236,
                                                                       1175, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3200, 0, 3, 6,
                                                                       2192, 698, 2210, 2948,
                                                                       1076, 2966, 236, 245,
                                                                       1202, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3254, 0, 3, 6,
                                                                       2210, 707, 2228, 2966,
                                                                       1085, 2984, 245, 254,
                                                                       1229, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3308, 0, 3, 6,
                                                                       2228, 716, 2246, 2984,
                                                                       1094, 3002, 254, 263,
                                                                       1256, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3362, 0, 3, 6,
                                                                       2246, 725, 2264, 3002,
                                                                       1103, 3020, 263, 272,
                                                                       1283, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 3416, 0, 3, 6,
                                                                       2282, 743, 2318, 3038,
                                                                       1121, 3092, 290, 308,
                                                                       1310, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 3524, 0, 3, 6,
                                                                       2318, 761, 2354, 3092,
                                                                       1148, 3146, 308, 326,
                                                                       1364, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 3632, 0, 3, 6,
                                                                       2354, 779, 2390, 3146,
                                                                       1175, 3200, 326, 344,
                                                                       1418, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 3740, 0, 3, 6,
                                                                       2390, 797, 2426, 3200,
                                                                       1202, 3254, 344, 362,
                                                                       1472, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 3848, 0, 3, 6,
                                                                       2426, 815, 2462, 3254,
                                                                       1229, 3308, 362, 380,
                                                                       1526, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 3956, 0, 3, 6,
                                                                       2462, 833, 2498, 3308,
                                                                       1256, 3362, 380, 398,
                                                                       1580, ncols, gamma, p,
                                                                       q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 4064, 0, 3, 6,
                                                                       2534, 869, 2594, 3038,
                                                                       3092, 3416, 1310, 3524,
                                                                       434, 464, 1634, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 4244, 0, 3, 6,
                                                                       2594, 899, 2654, 3092,
                                                                       3146, 3524, 1364, 3632,
                                                                       464, 494, 1724, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 4424, 0, 3, 6,
                                                                       2654, 929, 2714, 3146,
                                                                       3200, 3632, 1418, 3740,
                                                                       494, 524, 1814, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 4604, 0, 3, 6,
                                                                       2714, 959, 2774, 3200,
                                                                       3254, 3740, 1472, 3848,
                                                                       524, 554, 1904, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 4784, 0, 3, 6,
                                                                       2774, 989, 2834, 3254,
                                                                       3308, 3848, 1526, 3956,
                                                                       554, 584, 1994, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4964, 6, 644, 647,
                                                                       2096, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4974, 6, 647, 650,
                                                                       2102, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4984, 6, 650, 653,
                                                                       2108, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4994, 6, 653, 656,
                                                                       2114, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5004, 6, 656, 659,
                                                                       2120, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5014, 6, 659, 662,
                                                                       2126, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5024, 6, 662, 665,
                                                                       2132, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5034, 3, 6, 4964,
                                                                       2096, 4974, 2174, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5064, 3, 6, 4974,
                                                                       2102, 4984, 2192, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5094, 3, 6, 4984,
                                                                       2108, 4994, 2210, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5124, 3, 6, 4994,
                                                                       2114, 5004, 2228, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5154, 3, 6, 5004,
                                                                       2120, 5014, 2246, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5184, 3, 6, 5014,
                                                                       2126, 5024, 2264, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 5214, 3, 6, 5034,
                                                                       2174, 5064, 743, 761,
                                                                       2354, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 5274, 3, 6, 5064,
                                                                       2192, 5094, 761, 779,
                                                                       2390, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 5334, 3, 6, 5094,
                                                                       2210, 5124, 779, 797,
                                                                       2426, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 5394, 3, 6, 5124,
                                                                       2228, 5154, 797, 815,
                                                                       2462, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 5454, 3, 6, 5154,
                                                                       2246, 5184, 815, 833,
                                                                       2498, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 5514, 3, 6, 5214,
                                                                       2354, 5274, 869, 899,
                                                                       2654, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 5614, 3, 6, 5274,
                                                                       2390, 5334, 899, 929,
                                                                       2714, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 5714, 3, 6, 5334,
                                                                       2426, 5394, 929, 959,
                                                                       2774, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 5814, 3, 6, 5394,
                                                                       2462, 5454, 959, 989,
                                                                       2834, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5914, 0, 6, 4964,
                                                                       2096, 4974, 1049, 1058,
                                                                       2930, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5944, 0, 6, 4974,
                                                                       2102, 4984, 1058, 1067,
                                                                       2948, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5974, 0, 6, 4984,
                                                                       2108, 4994, 1067, 1076,
                                                                       2966, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6004, 0, 6, 4994,
                                                                       2114, 5004, 1076, 1085,
                                                                       2984, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6034, 0, 6, 5004,
                                                                       2120, 5014, 1085, 1094,
                                                                       3002, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6064, 0, 6, 5014,
                                                                       2126, 5024, 1094, 1103,
                                                                       3020, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 6094, 0, 3, 6,
                                                                       5034, 2174, 5064, 5914,
                                                                       2930, 5944, 1121, 1148,
                                                                       3146, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 6184, 0, 3, 6,
                                                                       5064, 2192, 5094, 5944,
                                                                       2948, 5974, 1148, 1175,
                                                                       3200, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 6274, 0, 3, 6,
                                                                       5094, 2210, 5124, 5974,
                                                                       2966, 6004, 1175, 1202,
                                                                       3254, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 6364, 0, 3, 6,
                                                                       5124, 2228, 5154, 6004,
                                                                       2984, 6034, 1202, 1229,
                                                                       3308, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 6454, 0, 3, 6,
                                                                       5154, 2246, 5184, 6034,
                                                                       3002, 6064, 1229, 1256,
                                                                       3362, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 6544, 0, 3, 6,
                                                                       5214, 2354, 5274, 6094,
                                                                       3146, 6184, 1310, 1364,
                                                                       3632, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 6724, 0, 3, 6,
                                                                       5274, 2390, 5334, 6184,
                                                                       3200, 6274, 1364, 1418,
                                                                       3740, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 6904, 0, 3, 6,
                                                                       5334, 2426, 5394, 6274,
                                                                       3254, 6364, 1418, 1472,
                                                                       3848, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 7084, 0, 3, 6,
                                                                       5394, 2462, 5454, 6364,
                                                                       3308, 6454, 1472, 1526,
                                                                       3956, ncols, gamma, p,
                                                                       q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 7264, 0, 3, 6,
                                                                       5514, 2654, 5614, 6094,
                                                                       6184, 6544, 3632, 6724,
                                                                       1634, 1724, 4424, ncols,
                                                                       gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 7564, 0, 3, 6,
                                                                       5614, 2714, 5714, 6184,
                                                                       6274, 6724, 3740, 6904,
                                                                       1724, 1814, 4604, ncols,
                                                                       gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 7864, 0, 3, 6,
                                                                       5714, 2774, 5814, 6274,
                                                                       6364, 6904, 3848, 7084,
                                                                       1814, 1904, 4784, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8164, 6, 2084,
                                                                       2090, 4964, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8179, 6, 2090,
                                                                       2096, 4974, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8194, 6, 2096,
                                                                       2102, 4984, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8209, 6, 2102,
                                                                       2108, 4994, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8224, 6, 2108,
                                                                       2114, 5004, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8239, 6, 2114,
                                                                       2120, 5014, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8254, 6, 2120,
                                                                       2126, 5024, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8269, 3, 6, 8164,
                                                                       4964, 8179, 2138, 2156,
                                                                       5034, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8314, 3, 6, 8179,
                                                                       4974, 8194, 2156, 2174,
                                                                       5064, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8359, 3, 6, 8194,
                                                                       4984, 8209, 2174, 2192,
                                                                       5094, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8404, 3, 6, 8209,
                                                                       4994, 8224, 2192, 2210,
                                                                       5124, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8449, 3, 6, 8224,
                                                                       5004, 8239, 2210, 2228,
                                                                       5154, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8494, 3, 6, 8239,
                                                                       5014, 8254, 2228, 2246,
                                                                       5184, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 8539, 3, 6, 8269,
                                                                       5034, 8314, 2282, 2318,
                                                                       5214, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 8629, 3, 6, 8314,
                                                                       5064, 8359, 2318, 2354,
                                                                       5274, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 8719, 3, 6, 8359,
                                                                       5094, 8404, 2354, 2390,
                                                                       5334, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 8809, 3, 6, 8404,
                                                                       5124, 8449, 2390, 2426,
                                                                       5394, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 8899, 3, 6, 8449,
                                                                       5154, 8494, 2426, 2462,
                                                                       5454, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 8989, 3, 6, 8539,
                                                                       5214, 8629, 2534, 2594,
                                                                       5514, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 9139, 3, 6, 8629,
                                                                       5274, 8719, 2594, 2654,
                                                                       5614, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 9289, 3, 6, 8719,
                                                                       5334, 8809, 2654, 2714,
                                                                       5714, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 9439, 3, 6, 8809,
                                                                       5394, 8899, 2714, 2774,
                                                                       5814, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 9589, 0, 6, 8164,
                                                                       4964, 8179, 2894, 2912,
                                                                       5914, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 9634, 0, 6, 8179,
                                                                       4974, 8194, 2912, 2930,
                                                                       5944, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 9679, 0, 6, 8194,
                                                                       4984, 8209, 2930, 2948,
                                                                       5974, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 9724, 0, 6, 8209,
                                                                       4994, 8224, 2948, 2966,
                                                                       6004, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 9769, 0, 6, 8224,
                                                                       5004, 8239, 2966, 2984,
                                                                       6034, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 9814, 0, 6, 8239,
                                                                       5014, 8254, 2984, 3002,
                                                                       6064, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 9859, 0, 3, 6,
                                                                       8269, 5034, 8314, 9589,
                                                                       5914, 9634, 3038, 3092,
                                                                       6094, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 9994, 0, 3, 6,
                                                                       8314, 5064, 8359, 9634,
                                                                       5944, 9679, 3092, 3146,
                                                                       6184, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 10129, 0, 3, 6,
                                                                       8359, 5094, 8404, 9679,
                                                                       5974, 9724, 3146, 3200,
                                                                       6274, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 10264, 0, 3, 6,
                                                                       8404, 5124, 8449, 9724,
                                                                       6004, 9769, 3200, 3254,
                                                                       6364, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 10399, 0, 3, 6,
                                                                       8449, 5154, 8494, 9769,
                                                                       6034, 9814, 3254, 3308,
                                                                       6454, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 10534, 0, 3, 6,
                                                                       8539, 5214, 8629, 9859,
                                                                       6094, 9994, 3416, 3524,
                                                                       6544, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 10804, 0, 3, 6,
                                                                       8629, 5274, 8719, 9994,
                                                                       6184, 10129, 3524, 3632,
                                                                       6724, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 11074, 0, 3, 6,
                                                                       8719, 5334, 8809, 10129,
                                                                       6274, 10264, 3632, 3740,
                                                                       6904, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 11344, 0, 3, 6,
                                                                       8809, 5394, 8899, 10264,
                                                                       6364, 10399, 3740, 3848,
                                                                       7084, ncols, gamma, p,
                                                                       q);

                    compute_prim_pfg_three_center_electron_repulsion_0(buffer, 11614, 0, 3, 6,
                                                                       8989, 5514, 9139, 9859,
                                                                       9994, 10534, 6544, 10804,
                                                                       4064, 4244, 7264, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfg_three_center_electron_repulsion_0(buffer, 12064, 0, 3, 6,
                                                                       9139, 5614, 9289, 9994,
                                                                       10129, 10804, 6724, 11074,
                                                                       4244, 4424, 7564, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfg_three_center_electron_repulsion_0(buffer, 12514, 0, 3, 6,
                                                                       9289, 5714, 9439, 10129,
                                                                       10264, 11074, 6904, 11344,
                                                                       4424, 4604, 7864, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 12964, 6, 4964,
                                                                       4974, 8194, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 12985, 6, 4974,
                                                                       4984, 8209, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 13006, 6, 4984,
                                                                       4994, 8224, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 13027, 6, 4994,
                                                                       5004, 8239, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 13048, 6, 5004,
                                                                       5014, 8254, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 13069, 3, 6,
                                                                       12964, 8194, 12985, 5034,
                                                                       5064, 8359, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 13132, 3, 6,
                                                                       12985, 8209, 13006, 5064,
                                                                       5094, 8404, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 13195, 3, 6,
                                                                       13006, 8224, 13027, 5094,
                                                                       5124, 8449, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 13258, 3, 6,
                                                                       13027, 8239, 13048, 5124,
                                                                       5154, 8494, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 13321, 3, 6,
                                                                       13069, 8359, 13132, 5214,
                                                                       5274, 8719, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 13447, 3, 6,
                                                                       13132, 8404, 13195, 5274,
                                                                       5334, 8809, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 13573, 3, 6,
                                                                       13195, 8449, 13258, 5334,
                                                                       5394, 8899, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 13699, 3, 6,
                                                                       13321, 8719, 13447, 5514,
                                                                       5614, 9289, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 13909, 3, 6,
                                                                       13447, 8809, 13573, 5614,
                                                                       5714, 9439, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 14119, 0, 6,
                                                                       12964, 8194, 12985, 5914,
                                                                       5944, 9679, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 14182, 0, 6,
                                                                       12985, 8209, 13006, 5944,
                                                                       5974, 9724, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 14245, 0, 6,
                                                                       13006, 8224, 13027, 5974,
                                                                       6004, 9769, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 14308, 0, 6,
                                                                       13027, 8239, 13048, 6004,
                                                                       6034, 9814, ncols, gamma,
                                                                       p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 14371, 0, 3, 6,
                                                                       13069, 8359, 13132, 14119,
                                                                       9679, 14182, 6094, 6184,
                                                                       10129, ncols, gamma, p,
                                                                       q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 14560, 0, 3, 6,
                                                                       13132, 8404, 13195, 14182,
                                                                       9724, 14245, 6184, 6274,
                                                                       10264, ncols, gamma, p,
                                                                       q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 14749, 0, 3, 6,
                                                                       13195, 8449, 13258, 14245,
                                                                       9769, 14308, 6274, 6364,
                                                                       10399, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdh_three_center_electron_repulsion_0(buffer, 14938, 0, 3, 6,
                                                                       13321, 8719, 13447, 14371,
                                                                       10129, 14560, 6544, 6724,
                                                                       11074, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdh_three_center_electron_repulsion_0(buffer, 15316, 0, 3, 6,
                                                                       13447, 8809, 13573, 14560,
                                                                       10264, 14749, 6724, 6904,
                                                                       11344, ncols, gamma, p,
                                                                       q);

                    compute_prim_pfh_three_center_electron_repulsion_0(buffer, 15694, 0, 3, 6,
                                                                       13699, 9289, 13909, 14371,
                                                                       14560, 14938, 11074,
                                                                       15316, 7264, 7564, 12514,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 16324, 6, 8164,
                                                                       8179, 12964, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 16352, 6, 8179,
                                                                       8194, 12985, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 16380, 6, 8194,
                                                                       8209, 13006, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 16408, 6, 8209,
                                                                       8224, 13027, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 16436, 6, 8224,
                                                                       8239, 13048, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 16464, 3, 6,
                                                                       16324, 12964, 16352, 8269,
                                                                       8314, 13069, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 16548, 3, 6,
                                                                       16352, 12985, 16380, 8314,
                                                                       8359, 13132, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 16632, 3, 6,
                                                                       16380, 13006, 16408, 8359,
                                                                       8404, 13195, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 16716, 3, 6,
                                                                       16408, 13027, 16436, 8404,
                                                                       8449, 13258, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 16800, 3, 6,
                                                                       16464, 13069, 16548, 8539,
                                                                       8629, 13321, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 16968, 3, 6,
                                                                       16548, 13132, 16632, 8629,
                                                                       8719, 13447, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 17136, 3, 6,
                                                                       16632, 13195, 16716, 8719,
                                                                       8809, 13573, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 17304, 3, 6,
                                                                       16800, 13321, 16968, 8989,
                                                                       9139, 13699, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 17584, 3, 6,
                                                                       16968, 13447, 17136, 9139,
                                                                       9289, 13909, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 17864, 0, 6,
                                                                       16324, 12964, 16352, 9589,
                                                                       9634, 14119, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 17948, 0, 6,
                                                                       16352, 12985, 16380, 9634,
                                                                       9679, 14182, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 18032, 0, 6,
                                                                       16380, 13006, 16408, 9679,
                                                                       9724, 14245, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 18116, 0, 6,
                                                                       16408, 13027, 16436, 9724,
                                                                       9769, 14308, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 18200, 0, 3, 6,
                                                                       16464, 13069, 16548,
                                                                       17864, 14119, 17948, 9859,
                                                                       9994, 14371, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 18452, 0, 3, 6,
                                                                       16548, 13132, 16632,
                                                                       17948, 14182, 18032, 9994,
                                                                       10129, 14560, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 18704, 0, 3, 6,
                                                                       16632, 13195, 16716,
                                                                       18032, 14245, 18116,
                                                                       10129, 10264, 14749,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdi_three_center_electron_repulsion_0(buffer, 18956, 0, 3, 6,
                                                                       16800, 13321, 16968,
                                                                       18200, 14371, 18452,
                                                                       10534, 10804, 14938,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdi_three_center_electron_repulsion_0(buffer, 19460, 0, 3, 6,
                                                                       16968, 13447, 17136,
                                                                       18452, 14560, 18704,
                                                                       10804, 11074, 15316,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfi_three_center_electron_repulsion_0(buffer, 19964, 0, 3, 6,
                                                                       17304, 13699, 17584,
                                                                       18200, 18452, 18956,
                                                                       14938, 19460, 11614,
                                                                       12064, 15694, ncols,
                                                                       gamma, p, q);

                    simdgeo::geom_s_x(buffer, 20804, 19964, 1, 280, ncols, alpha);

                    simdgeo::geom_s_y(buffer, 21084, 19964, 1, 280, ncols, alpha);

                    simdgeo::geom_s_z(buffer, 21364, 19964, 1, 280, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 21644, 20804, 840, ncols);
                }
            }
        }

        simdtrf::transform_i_inner(buffer, 22484, 21644, 10, 1, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 22484, 13, nmax);

        simdtrf::transform_i_inner(buffer, 22484, 21924, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 91 * nvalues + n * npairs, nvalues, buffer, 22484,
                                   13, nmax);

        simdtrf::transform_i_inner(buffer, 22484, 22204, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 182 * nvalues + n * npairs, nvalues, buffer, 22484,
                                   13, nmax);
    }

    for (size_t m = 0; m < 273; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
