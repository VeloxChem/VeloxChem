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


#include "SimdThreeCenterElectronRepulsionGeom100RecSGF.hpp"

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
compute_geom_100_sgf_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_geom_100_sgf_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 8417, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 189 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 8417, 7862, 450, dimensions);

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

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 9, 6, {1, 2, 3, 4,
                                                        5, 6, 7, 8}, ncols, fj, i * nprim_b + j,
                                                        fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 18, 3, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 21, 3, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 24, 3, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 27, 3, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 30, 3, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 33, 3, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 36, 3, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 39, 3, 6, 10, 11,
                                                                       18, 21, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 45, 3, 6, 11, 12,
                                                                       21, 24, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 51, 3, 6, 12, 13,
                                                                       24, 27, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 57, 3, 6, 13, 14,
                                                                       27, 30, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 63, 3, 6, 14, 15,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 69, 3, 6, 15, 16,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 75, 3, 6, 18, 21,
                                                                       39, 45, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 85, 3, 6, 21, 24,
                                                                       45, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 95, 3, 6, 24, 27,
                                                                       51, 57, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 105, 3, 6, 27, 30,
                                                                       57, 63, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 115, 3, 6, 30, 33,
                                                                       63, 69, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 125, 3, 6, 39, 45,
                                                                       75, 85, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 140, 3, 6, 45, 51,
                                                                       85, 95, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 155, 3, 6, 51, 57,
                                                                       95, 105, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 170, 3, 6, 57, 63,
                                                                       105, 115, ncols, gamma, p,
                                                                       q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 185, 0, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 188, 0, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 191, 0, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 194, 0, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 197, 0, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 200, 0, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 203, 0, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 206, 0, 6, 10, 11,
                                                                       18, 21, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 215, 0, 6, 11, 12,
                                                                       21, 24, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 224, 0, 6, 12, 13,
                                                                       24, 27, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 233, 0, 6, 13, 14,
                                                                       27, 30, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 242, 0, 6, 14, 15,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 251, 0, 6, 15, 16,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 260, 0, 3, 6, 18,
                                                                       21, 39, 45, 206, 215,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 278, 0, 3, 6, 21,
                                                                       24, 45, 51, 215, 224,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 296, 0, 3, 6, 24,
                                                                       27, 51, 57, 224, 233,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 314, 0, 3, 6, 27,
                                                                       30, 57, 63, 233, 242,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 332, 0, 3, 6, 30,
                                                                       33, 63, 69, 242, 251,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 350, 0, 3, 6, 39,
                                                                       45, 75, 85, 260, 278,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 380, 0, 3, 6, 45,
                                                                       51, 85, 95, 278, 296,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 410, 0, 3, 6, 51,
                                                                       57, 95, 105, 296, 314,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 440, 0, 3, 6, 57,
                                                                       63, 105, 115, 314, 332,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 470, 0, 3, 6, 75,
                                                                       85, 125, 140, 350, 380,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 515, 0, 3, 6, 85,
                                                                       95, 140, 155, 380, 410,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 560, 0, 3, 6, 95,
                                                                       105, 155, 170, 410, 440,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 605, 6, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 608, 6, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 611, 6, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 614, 6, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 617, 6, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 620, 6, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 623, 6, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 626, 6, 17, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 629, 6, 12, 24,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 638, 6, 13, 27,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 647, 6, 14, 30,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 656, 6, 15, 33,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 665, 6, 16, 36,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 674, 6, 18, 39,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 692, 6, 21, 45,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 710, 6, 24, 51,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 728, 6, 27, 57,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 746, 6, 30, 63,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 764, 6, 33, 69,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 782, 6, 39, 75,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 812, 6, 45, 85,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 842, 6, 51, 95,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 872, 6, 57, 105,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 902, 6, 63, 115,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 932, 6, 75, 125,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 977, 6, 85, 140,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1022, 6, 95, 155,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1067, 6, 105, 170,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1112, 6, 10, 185,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1121, 6, 11, 188,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1130, 6, 12, 191,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1139, 6, 13, 194,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1148, 6, 14, 197,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1157, 6, 15, 200,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1166, 6, 16, 203,
                                                                       ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1175, 6, 18, 185,
                                                                       206, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1202, 6, 21, 188,
                                                                       215, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1229, 6, 24, 191,
                                                                       224, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1256, 6, 27, 194,
                                                                       233, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1283, 6, 30, 197,
                                                                       242, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1310, 6, 33, 200,
                                                                       251, ncols, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1337, 3, 6, 39,
                                                                       1175, 206, 1202, 260,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1391, 3, 6, 45,
                                                                       1202, 215, 1229, 278,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1445, 3, 6, 51,
                                                                       1229, 224, 1256, 296,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1499, 3, 6, 57,
                                                                       1256, 233, 1283, 314,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1553, 3, 6, 63,
                                                                       1283, 242, 1310, 332,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 1607, 3, 6, 75,
                                                                       1337, 260, 1391, 350,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 1697, 3, 6, 85,
                                                                       1391, 278, 1445, 380,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 1787, 3, 6, 95,
                                                                       1445, 296, 1499, 410,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 1877, 3, 6, 105,
                                                                       1499, 314, 1553, 440,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 1967, 3, 6, 125,
                                                                       1607, 350, 1697, 470,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 2102, 3, 6, 140,
                                                                       1697, 380, 1787, 515,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 2237, 3, 6, 155,
                                                                       1787, 410, 1877, 560,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2372, 6, 10, 11,
                                                                       611, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2378, 6, 11, 12,
                                                                       614, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2384, 6, 12, 13,
                                                                       617, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2390, 6, 13, 14,
                                                                       620, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2396, 6, 14, 15,
                                                                       623, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2402, 6, 15, 16,
                                                                       626, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2408, 3, 6, 2372,
                                                                       611, 2378, 629, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2426, 3, 6, 2378,
                                                                       614, 2384, 638, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2444, 3, 6, 2384,
                                                                       617, 2390, 647, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2462, 3, 6, 2390,
                                                                       620, 2396, 656, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2480, 3, 6, 2396,
                                                                       623, 2402, 665, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2498, 3, 6, 2408,
                                                                       629, 2426, 39, 45, 710,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2534, 3, 6, 2426,
                                                                       638, 2444, 45, 51, 728,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2570, 3, 6, 2444,
                                                                       647, 2462, 51, 57, 746,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2606, 3, 6, 2462,
                                                                       656, 2480, 57, 63, 764,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2642, 3, 6, 2498,
                                                                       710, 2534, 75, 85, 842,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2702, 3, 6, 2534,
                                                                       728, 2570, 85, 95, 872,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2762, 3, 6, 2570,
                                                                       746, 2606, 95, 105, 902,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 2822, 3, 6, 2642,
                                                                       842, 2702, 125, 140, 1022,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 2912, 3, 6, 2702,
                                                                       872, 2762, 140, 155, 1067,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3002, 0, 6, 2372,
                                                                       611, 2378, 1130, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3020, 0, 6, 2378,
                                                                       614, 2384, 1139, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3038, 0, 6, 2384,
                                                                       617, 2390, 1148, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3056, 0, 6, 2390,
                                                                       620, 2396, 1157, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3074, 0, 6, 2396,
                                                                       623, 2402, 1166, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3092, 0, 3, 6,
                                                                       2408, 629, 2426, 3002,
                                                                       1130, 3020, 206, 215,
                                                                       1229, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3146, 0, 3, 6,
                                                                       2426, 638, 2444, 3020,
                                                                       1139, 3038, 215, 224,
                                                                       1256, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3200, 0, 3, 6,
                                                                       2444, 647, 2462, 3038,
                                                                       1148, 3056, 224, 233,
                                                                       1283, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3254, 0, 3, 6,
                                                                       2462, 656, 2480, 3056,
                                                                       1157, 3074, 233, 242,
                                                                       1310, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 3308, 0, 3, 6,
                                                                       2498, 710, 2534, 3092,
                                                                       1229, 3146, 260, 278,
                                                                       1445, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 3416, 0, 3, 6,
                                                                       2534, 728, 2570, 3146,
                                                                       1256, 3200, 278, 296,
                                                                       1499, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 3524, 0, 3, 6,
                                                                       2570, 746, 2606, 3200,
                                                                       1283, 3254, 296, 314,
                                                                       1553, ncols, gamma, p,
                                                                       q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 3632, 0, 3, 6,
                                                                       2642, 842, 2702, 3092,
                                                                       3146, 3308, 1445, 3416,
                                                                       350, 380, 1787, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 3812, 0, 3, 6,
                                                                       2702, 872, 2762, 3146,
                                                                       3200, 3416, 1499, 3524,
                                                                       380, 410, 1877, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgd_three_center_electron_repulsion_0(buffer, 3992, 0, 3, 6,
                                                                       2822, 1022, 2912, 3308,
                                                                       3416, 3632, 1787, 3812,
                                                                       470, 515, 2237, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4262, 6, 605, 608,
                                                                       2372, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4272, 6, 608, 611,
                                                                       2378, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4282, 6, 611, 614,
                                                                       2384, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4292, 6, 614, 617,
                                                                       2390, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4302, 6, 617, 620,
                                                                       2396, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4312, 6, 620, 623,
                                                                       2402, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4322, 3, 6, 4262,
                                                                       2372, 4272, 2408, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4352, 3, 6, 4272,
                                                                       2378, 4282, 2426, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4382, 3, 6, 4282,
                                                                       2384, 4292, 2444, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4412, 3, 6, 4292,
                                                                       2390, 4302, 2462, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4442, 3, 6, 4302,
                                                                       2396, 4312, 2480, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4472, 3, 6, 4322,
                                                                       2408, 4352, 674, 692,
                                                                       2498, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4532, 3, 6, 4352,
                                                                       2426, 4382, 692, 710,
                                                                       2534, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4592, 3, 6, 4382,
                                                                       2444, 4412, 710, 728,
                                                                       2570, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4652, 3, 6, 4412,
                                                                       2462, 4442, 728, 746,
                                                                       2606, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 4712, 3, 6, 4472,
                                                                       2498, 4532, 782, 812,
                                                                       2642, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 4812, 3, 6, 4532,
                                                                       2534, 4592, 812, 842,
                                                                       2702, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 4912, 3, 6, 4592,
                                                                       2570, 4652, 842, 872,
                                                                       2762, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 5012, 3, 6, 4712,
                                                                       2642, 4812, 932, 977,
                                                                       2822, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 5162, 3, 6, 4812,
                                                                       2702, 4912, 977, 1022,
                                                                       2912, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5312, 0, 6, 4262,
                                                                       2372, 4272, 1112, 1121,
                                                                       3002, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5342, 0, 6, 4272,
                                                                       2378, 4282, 1121, 1130,
                                                                       3020, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5372, 0, 6, 4282,
                                                                       2384, 4292, 1130, 1139,
                                                                       3038, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5402, 0, 6, 4292,
                                                                       2390, 4302, 1139, 1148,
                                                                       3056, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5432, 0, 6, 4302,
                                                                       2396, 4312, 1148, 1157,
                                                                       3074, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 5462, 0, 3, 6,
                                                                       4322, 2408, 4352, 5312,
                                                                       3002, 5342, 1175, 1202,
                                                                       3092, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 5552, 0, 3, 6,
                                                                       4352, 2426, 4382, 5342,
                                                                       3020, 5372, 1202, 1229,
                                                                       3146, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 5642, 0, 3, 6,
                                                                       4382, 2444, 4412, 5372,
                                                                       3038, 5402, 1229, 1256,
                                                                       3200, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 5732, 0, 3, 6,
                                                                       4412, 2462, 4442, 5402,
                                                                       3056, 5432, 1256, 1283,
                                                                       3254, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 5822, 0, 3, 6,
                                                                       4472, 2498, 4532, 5462,
                                                                       3092, 5552, 1337, 1391,
                                                                       3308, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 6002, 0, 3, 6,
                                                                       4532, 2534, 4592, 5552,
                                                                       3146, 5642, 1391, 1445,
                                                                       3416, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 6182, 0, 3, 6,
                                                                       4592, 2570, 4652, 5642,
                                                                       3200, 5732, 1445, 1499,
                                                                       3524, ncols, gamma, p,
                                                                       q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 6362, 0, 3, 6,
                                                                       4712, 2642, 4812, 5462,
                                                                       5552, 5822, 3308, 6002,
                                                                       1607, 1697, 3632, ncols,
                                                                       gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 6662, 0, 3, 6,
                                                                       4812, 2702, 4912, 5552,
                                                                       5642, 6002, 3416, 6182,
                                                                       1697, 1787, 3812, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgf_three_center_electron_repulsion_0(buffer, 6962, 0, 3, 6,
                                                                       5012, 2822, 5162, 5822,
                                                                       6002, 6362, 3632, 6662,
                                                                       1967, 2102, 3992, ncols,
                                                                       gamma, p, q);

                    simdgeo::geom_s_x(buffer, 7412, 6962, 1, 150, ncols, alpha);

                    simdgeo::geom_s_y(buffer, 7562, 6962, 1, 150, ncols, alpha);

                    simdgeo::geom_s_z(buffer, 7712, 6962, 1, 150, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 7862, 7412, 450, ncols);
                }
            }
        }

        simdtrf::transform_f_inner(buffer, 8312, 7862, 15, 1, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 8312, 7, nmax);

        simdtrf::transform_f_inner(buffer, 8312, 8012, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 63 * nvalues + n * npairs, nvalues, buffer, 8312, 7,
                                   nmax);

        simdtrf::transform_f_inner(buffer, 8312, 8162, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 126 * nvalues + n * npairs, nvalues, buffer, 8312, 7,
                                   nmax);
    }

    for (size_t m = 0; m < 189; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
