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


#include "SimdThreeCenterElectronRepulsionGeom100RecSDG.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformG.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_geom_100_sdg_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_geom_100_sdg_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 3928, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 135 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 3928, 3604, 270, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 9, 6, 7, ncols,
                                                             fj, i * nprim_b + j, fq);

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

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 75, 0, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 78, 0, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 81, 0, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 84, 0, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 87, 0, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 90, 0, 6, 10, 11,
                                                                       18, 21, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 99, 0, 6, 11, 12,
                                                                       21, 24, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 108, 0, 6, 12, 13,
                                                                       24, 27, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 117, 0, 6, 13, 14,
                                                                       27, 30, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 126, 0, 6, 14, 15,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 135, 0, 6, 15, 16,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 144, 0, 3, 6, 18,
                                                                       21, 39, 45, 90, 99, ncols,
                                                                       gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 162, 0, 3, 6, 21,
                                                                       24, 45, 51, 99, 108,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 180, 0, 3, 6, 24,
                                                                       27, 51, 57, 108, 117,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 198, 0, 3, 6, 27,
                                                                       30, 57, 63, 117, 126,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 216, 0, 3, 6, 30,
                                                                       33, 63, 69, 126, 135,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 234, 6, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 237, 6, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 240, 6, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 243, 6, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 246, 6, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 249, 6, 17, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 252, 6, 12, 24,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 261, 6, 13, 27,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 270, 6, 14, 30,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 279, 6, 15, 33,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 288, 6, 16, 36,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 297, 6, 24, 51,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 315, 6, 27, 57,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 333, 6, 30, 63,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 351, 6, 33, 69,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 369, 6, 12, 75,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 378, 6, 13, 78,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 387, 6, 14, 81,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 396, 6, 15, 84,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 405, 6, 16, 87,
                                                                       ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 414, 6, 24, 75,
                                                                       108, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 441, 6, 27, 78,
                                                                       117, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 468, 6, 30, 81,
                                                                       126, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 495, 6, 33, 84,
                                                                       135, ncols, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 522, 3, 6, 51,
                                                                       414, 108, 441, 180, ncols,
                                                                       gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 576, 3, 6, 57,
                                                                       441, 117, 468, 198, ncols,
                                                                       gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 630, 3, 6, 63,
                                                                       468, 126, 495, 216, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 684, 6, 10, 11,
                                                                       234, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 690, 6, 11, 12,
                                                                       237, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 696, 6, 12, 13,
                                                                       240, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 702, 6, 13, 14,
                                                                       243, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 708, 6, 14, 15,
                                                                       246, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 714, 6, 15, 16,
                                                                       249, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 720, 3, 6, 684,
                                                                       234, 690, 252, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 738, 3, 6, 690,
                                                                       237, 696, 261, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 756, 3, 6, 696,
                                                                       240, 702, 270, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 774, 3, 6, 702,
                                                                       243, 708, 279, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 792, 3, 6, 708,
                                                                       246, 714, 288, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 810, 3, 6, 720,
                                                                       252, 738, 39, 45, 297,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 846, 3, 6, 738,
                                                                       261, 756, 45, 51, 315,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 882, 3, 6, 756,
                                                                       270, 774, 51, 57, 333,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 918, 3, 6, 774,
                                                                       279, 792, 57, 63, 351,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 954, 0, 6, 684,
                                                                       234, 690, 369, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 972, 0, 6, 690,
                                                                       237, 696, 378, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 990, 0, 6, 696,
                                                                       240, 702, 387, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1008, 0, 6, 702,
                                                                       243, 708, 396, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1026, 0, 6, 708,
                                                                       246, 714, 405, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 1044, 0, 3, 6,
                                                                       720, 252, 738, 954, 369,
                                                                       972, 90, 99, 414, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 1098, 0, 3, 6,
                                                                       738, 261, 756, 972, 378,
                                                                       990, 99, 108, 441, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 1152, 0, 3, 6,
                                                                       756, 270, 774, 990, 387,
                                                                       1008, 108, 117, 468,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 1206, 0, 3, 6,
                                                                       774, 279, 792, 1008, 396,
                                                                       1026, 117, 126, 495,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 1260, 0, 3, 6,
                                                                       810, 297, 846, 1044, 414,
                                                                       1098, 144, 162, 522,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 1368, 0, 3, 6,
                                                                       846, 315, 882, 1098, 441,
                                                                       1152, 162, 180, 576,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 1476, 0, 3, 6,
                                                                       882, 333, 918, 1152, 468,
                                                                       1206, 180, 198, 630,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1584, 6, 234, 237,
                                                                       696, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1594, 6, 237, 240,
                                                                       702, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1604, 6, 240, 243,
                                                                       708, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1614, 6, 243, 246,
                                                                       714, ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 1624, 3, 6, 1584,
                                                                       696, 1594, 756, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 1654, 3, 6, 1594,
                                                                       702, 1604, 774, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 1684, 3, 6, 1604,
                                                                       708, 1614, 792, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 1714, 3, 6, 1624,
                                                                       756, 1654, 297, 315, 882,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 1774, 3, 6, 1654,
                                                                       774, 1684, 315, 333, 918,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 1834, 0, 6, 1584,
                                                                       696, 1594, 369, 378, 990,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 1864, 0, 6, 1594,
                                                                       702, 1604, 378, 387, 1008,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 1894, 0, 6, 1604,
                                                                       708, 1614, 387, 396, 1026,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 1924, 0, 3, 6,
                                                                       1624, 756, 1654, 1834,
                                                                       990, 1864, 414, 441, 1152,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 2014, 0, 3, 6,
                                                                       1654, 774, 1684, 1864,
                                                                       1008, 1894, 441, 468,
                                                                       1206, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 2104, 0, 3, 6,
                                                                       1714, 882, 1774, 1924,
                                                                       1152, 2014, 522, 576,
                                                                       1476, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 2284, 6, 684, 690,
                                                                       1584, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 2299, 6, 690, 696,
                                                                       1594, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 2314, 6, 696, 702,
                                                                       1604, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 2329, 6, 702, 708,
                                                                       1614, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 2344, 3, 6, 2284,
                                                                       1584, 2299, 720, 738,
                                                                       1624, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 2389, 3, 6, 2299,
                                                                       1594, 2314, 738, 756,
                                                                       1654, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 2434, 3, 6, 2314,
                                                                       1604, 2329, 756, 774,
                                                                       1684, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 2479, 3, 6, 2344,
                                                                       1624, 2389, 810, 846,
                                                                       1714, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 2569, 3, 6, 2389,
                                                                       1654, 2434, 846, 882,
                                                                       1774, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 2659, 0, 6, 2284,
                                                                       1584, 2299, 954, 972,
                                                                       1834, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 2704, 0, 6, 2299,
                                                                       1594, 2314, 972, 990,
                                                                       1864, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 2749, 0, 6, 2314,
                                                                       1604, 2329, 990, 1008,
                                                                       1894, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 2794, 0, 3, 6,
                                                                       2344, 1624, 2389, 2659,
                                                                       1834, 2704, 1044, 1098,
                                                                       1924, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 2929, 0, 3, 6,
                                                                       2389, 1654, 2434, 2704,
                                                                       1864, 2749, 1098, 1152,
                                                                       2014, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 3064, 0, 3, 6,
                                                                       2479, 1714, 2569, 2794,
                                                                       1924, 2929, 1260, 1368,
                                                                       2104, ncols, gamma, p,
                                                                       q);

                    simdgeo::geom_s_x(buffer, 3334, 3064, 1, 90, ncols, alpha);

                    simdgeo::geom_s_y(buffer, 3424, 3064, 1, 90, ncols, alpha);

                    simdgeo::geom_s_z(buffer, 3514, 3064, 1, 90, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 3604, 3334, 270, ncols);
                }
            }
        }

        simdtrf::transform_g_inner(buffer, 3874, 3604, 6, 1, nmax);

        simdtrf::transform_d_outer(values + n * npairs, nvalues, buffer, 3874, 9, nmax);

        simdtrf::transform_g_inner(buffer, 3874, 3694, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 45 * nvalues + n * npairs, nvalues, buffer, 3874, 9,
                                   nmax);

        simdtrf::transform_g_inner(buffer, 3874, 3784, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 90 * nvalues + n * npairs, nvalues, buffer, 3874, 9,
                                   nmax);
    }

    for (size_t m = 0; m < 135; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
