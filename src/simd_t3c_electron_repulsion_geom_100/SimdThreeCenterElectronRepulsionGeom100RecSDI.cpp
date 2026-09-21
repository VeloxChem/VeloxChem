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


#include "SimdThreeCenterElectronRepulsionGeom100RecSDI.hpp"

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
compute_geom_100_sdi_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_geom_100_sdi_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 10650, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 195 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 10650, 10068, 504, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 9, 6, 9, ncols,
                                                             fj, i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 20, 3, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 23, 3, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 26, 3, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 29, 3, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 32, 3, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 35, 3, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 38, 3, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 41, 3, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 44, 3, 6, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 47, 3, 6, 10, 11,
                                                                       20, 23, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 53, 3, 6, 11, 12,
                                                                       23, 26, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 59, 3, 6, 12, 13,
                                                                       26, 29, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 65, 3, 6, 13, 14,
                                                                       29, 32, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 71, 3, 6, 14, 15,
                                                                       32, 35, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 77, 3, 6, 15, 16,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 83, 3, 6, 16, 17,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 89, 3, 6, 17, 18,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 95, 0, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 98, 0, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 101, 0, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 104, 0, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 107, 0, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 110, 0, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 113, 0, 6, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 116, 0, 6, 10, 11,
                                                                       20, 23, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 125, 0, 6, 11, 12,
                                                                       23, 26, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 134, 0, 6, 12, 13,
                                                                       26, 29, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 143, 0, 6, 13, 14,
                                                                       29, 32, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 152, 0, 6, 14, 15,
                                                                       32, 35, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 161, 0, 6, 15, 16,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 170, 0, 6, 16, 17,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 179, 0, 6, 17, 18,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 188, 0, 3, 6, 20,
                                                                       23, 47, 53, 116, 125,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 206, 0, 3, 6, 23,
                                                                       26, 53, 59, 125, 134,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 224, 0, 3, 6, 26,
                                                                       29, 59, 65, 134, 143,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 242, 0, 3, 6, 29,
                                                                       32, 65, 71, 143, 152,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 260, 0, 3, 6, 32,
                                                                       35, 71, 77, 152, 161,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 278, 0, 3, 6, 35,
                                                                       38, 77, 83, 161, 170,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 296, 0, 3, 6, 38,
                                                                       41, 83, 89, 170, 179,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 314, 6, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 317, 6, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 320, 6, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 323, 6, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 326, 6, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 329, 6, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 332, 6, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 335, 6, 19, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 338, 6, 12, 26,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 347, 6, 13, 29,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 356, 6, 14, 32,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 365, 6, 15, 35,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 374, 6, 16, 38,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 383, 6, 17, 41,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 392, 6, 18, 44,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 401, 6, 26, 59,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 419, 6, 29, 65,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 437, 6, 32, 71,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 455, 6, 35, 77,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 473, 6, 38, 83,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 491, 6, 41, 89,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 509, 6, 12, 95,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 518, 6, 13, 98,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 527, 6, 14, 101,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 536, 6, 15, 104,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 545, 6, 16, 107,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 554, 6, 17, 110,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 563, 6, 18, 113,
                                                                       ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 572, 6, 26, 95,
                                                                       134, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 599, 6, 29, 98,
                                                                       143, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 626, 6, 32, 101,
                                                                       152, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 653, 6, 35, 104,
                                                                       161, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 680, 6, 38, 107,
                                                                       170, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 707, 6, 41, 110,
                                                                       179, ncols, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 734, 3, 6, 59,
                                                                       572, 134, 599, 224, ncols,
                                                                       gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 788, 3, 6, 65,
                                                                       599, 143, 626, 242, ncols,
                                                                       gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 842, 3, 6, 71,
                                                                       626, 152, 653, 260, ncols,
                                                                       gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 896, 3, 6, 77,
                                                                       653, 161, 680, 278, ncols,
                                                                       gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 950, 3, 6, 83,
                                                                       680, 170, 707, 296, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1004, 6, 10, 11,
                                                                       314, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1010, 6, 11, 12,
                                                                       317, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1016, 6, 12, 13,
                                                                       320, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1022, 6, 13, 14,
                                                                       323, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1028, 6, 14, 15,
                                                                       326, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1034, 6, 15, 16,
                                                                       329, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1040, 6, 16, 17,
                                                                       332, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1046, 6, 17, 18,
                                                                       335, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1052, 3, 6, 1004,
                                                                       314, 1010, 338, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1070, 3, 6, 1010,
                                                                       317, 1016, 347, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1088, 3, 6, 1016,
                                                                       320, 1022, 356, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1106, 3, 6, 1022,
                                                                       323, 1028, 365, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1124, 3, 6, 1028,
                                                                       326, 1034, 374, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1142, 3, 6, 1034,
                                                                       329, 1040, 383, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1160, 3, 6, 1040,
                                                                       332, 1046, 392, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1178, 3, 6, 1052,
                                                                       338, 1070, 47, 53, 401,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1214, 3, 6, 1070,
                                                                       347, 1088, 53, 59, 419,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1250, 3, 6, 1088,
                                                                       356, 1106, 59, 65, 437,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1286, 3, 6, 1106,
                                                                       365, 1124, 65, 71, 455,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1322, 3, 6, 1124,
                                                                       374, 1142, 71, 77, 473,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1358, 3, 6, 1142,
                                                                       383, 1160, 77, 83, 491,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1394, 0, 6, 1004,
                                                                       314, 1010, 509, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1412, 0, 6, 1010,
                                                                       317, 1016, 518, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1430, 0, 6, 1016,
                                                                       320, 1022, 527, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1448, 0, 6, 1022,
                                                                       323, 1028, 536, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1466, 0, 6, 1028,
                                                                       326, 1034, 545, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1484, 0, 6, 1034,
                                                                       329, 1040, 554, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1502, 0, 6, 1040,
                                                                       332, 1046, 563, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 1520, 0, 3, 6,
                                                                       1052, 338, 1070, 1394,
                                                                       509, 1412, 116, 125, 572,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 1574, 0, 3, 6,
                                                                       1070, 347, 1088, 1412,
                                                                       518, 1430, 125, 134, 599,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 1628, 0, 3, 6,
                                                                       1088, 356, 1106, 1430,
                                                                       527, 1448, 134, 143, 626,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 1682, 0, 3, 6,
                                                                       1106, 365, 1124, 1448,
                                                                       536, 1466, 143, 152, 653,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 1736, 0, 3, 6,
                                                                       1124, 374, 1142, 1466,
                                                                       545, 1484, 152, 161, 680,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 1790, 0, 3, 6,
                                                                       1142, 383, 1160, 1484,
                                                                       554, 1502, 161, 170, 707,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 1844, 0, 3, 6,
                                                                       1178, 401, 1214, 1520,
                                                                       572, 1574, 188, 206, 734,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 1952, 0, 3, 6,
                                                                       1214, 419, 1250, 1574,
                                                                       599, 1628, 206, 224, 788,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 2060, 0, 3, 6,
                                                                       1250, 437, 1286, 1628,
                                                                       626, 1682, 224, 242, 842,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 2168, 0, 3, 6,
                                                                       1286, 455, 1322, 1682,
                                                                       653, 1736, 242, 260, 896,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 2276, 0, 3, 6,
                                                                       1322, 473, 1358, 1736,
                                                                       680, 1790, 260, 278, 950,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2384, 6, 314, 317,
                                                                       1016, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2394, 6, 317, 320,
                                                                       1022, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2404, 6, 320, 323,
                                                                       1028, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2414, 6, 323, 326,
                                                                       1034, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2424, 6, 326, 329,
                                                                       1040, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2434, 6, 329, 332,
                                                                       1046, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 2444, 3, 6, 2384,
                                                                       1016, 2394, 1088, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 2474, 3, 6, 2394,
                                                                       1022, 2404, 1106, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 2504, 3, 6, 2404,
                                                                       1028, 2414, 1124, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 2534, 3, 6, 2414,
                                                                       1034, 2424, 1142, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 2564, 3, 6, 2424,
                                                                       1040, 2434, 1160, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 2594, 3, 6, 2444,
                                                                       1088, 2474, 401, 419,
                                                                       1250, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 2654, 3, 6, 2474,
                                                                       1106, 2504, 419, 437,
                                                                       1286, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 2714, 3, 6, 2504,
                                                                       1124, 2534, 437, 455,
                                                                       1322, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 2774, 3, 6, 2534,
                                                                       1142, 2564, 455, 473,
                                                                       1358, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2834, 0, 6, 2384,
                                                                       1016, 2394, 509, 518,
                                                                       1430, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2864, 0, 6, 2394,
                                                                       1022, 2404, 518, 527,
                                                                       1448, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2894, 0, 6, 2404,
                                                                       1028, 2414, 527, 536,
                                                                       1466, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2924, 0, 6, 2414,
                                                                       1034, 2424, 536, 545,
                                                                       1484, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2954, 0, 6, 2424,
                                                                       1040, 2434, 545, 554,
                                                                       1502, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 2984, 0, 3, 6,
                                                                       2444, 1088, 2474, 2834,
                                                                       1430, 2864, 572, 599,
                                                                       1628, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 3074, 0, 3, 6,
                                                                       2474, 1106, 2504, 2864,
                                                                       1448, 2894, 599, 626,
                                                                       1682, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 3164, 0, 3, 6,
                                                                       2504, 1124, 2534, 2894,
                                                                       1466, 2924, 626, 653,
                                                                       1736, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 3254, 0, 3, 6,
                                                                       2534, 1142, 2564, 2924,
                                                                       1484, 2954, 653, 680,
                                                                       1790, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 3344, 0, 3, 6,
                                                                       2594, 1250, 2654, 2984,
                                                                       1628, 3074, 734, 788,
                                                                       2060, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 3524, 0, 3, 6,
                                                                       2654, 1286, 2714, 3074,
                                                                       1682, 3164, 788, 842,
                                                                       2168, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 3704, 0, 3, 6,
                                                                       2714, 1322, 2774, 3164,
                                                                       1736, 3254, 842, 896,
                                                                       2276, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 3884, 6, 1004,
                                                                       1010, 2384, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 3899, 6, 1010,
                                                                       1016, 2394, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 3914, 6, 1016,
                                                                       1022, 2404, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 3929, 6, 1022,
                                                                       1028, 2414, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 3944, 6, 1028,
                                                                       1034, 2424, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 3959, 6, 1034,
                                                                       1040, 2434, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 3974, 3, 6, 3884,
                                                                       2384, 3899, 1052, 1070,
                                                                       2444, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 4019, 3, 6, 3899,
                                                                       2394, 3914, 1070, 1088,
                                                                       2474, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 4064, 3, 6, 3914,
                                                                       2404, 3929, 1088, 1106,
                                                                       2504, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 4109, 3, 6, 3929,
                                                                       2414, 3944, 1106, 1124,
                                                                       2534, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 4154, 3, 6, 3944,
                                                                       2424, 3959, 1124, 1142,
                                                                       2564, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 4199, 3, 6, 3974,
                                                                       2444, 4019, 1178, 1214,
                                                                       2594, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 4289, 3, 6, 4019,
                                                                       2474, 4064, 1214, 1250,
                                                                       2654, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 4379, 3, 6, 4064,
                                                                       2504, 4109, 1250, 1286,
                                                                       2714, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 4469, 3, 6, 4109,
                                                                       2534, 4154, 1286, 1322,
                                                                       2774, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 4559, 0, 6, 3884,
                                                                       2384, 3899, 1394, 1412,
                                                                       2834, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 4604, 0, 6, 3899,
                                                                       2394, 3914, 1412, 1430,
                                                                       2864, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 4649, 0, 6, 3914,
                                                                       2404, 3929, 1430, 1448,
                                                                       2894, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 4694, 0, 6, 3929,
                                                                       2414, 3944, 1448, 1466,
                                                                       2924, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 4739, 0, 6, 3944,
                                                                       2424, 3959, 1466, 1484,
                                                                       2954, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 4784, 0, 3, 6,
                                                                       3974, 2444, 4019, 4559,
                                                                       2834, 4604, 1520, 1574,
                                                                       2984, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 4919, 0, 3, 6,
                                                                       4019, 2474, 4064, 4604,
                                                                       2864, 4649, 1574, 1628,
                                                                       3074, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 5054, 0, 3, 6,
                                                                       4064, 2504, 4109, 4649,
                                                                       2894, 4694, 1628, 1682,
                                                                       3164, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 5189, 0, 3, 6,
                                                                       4109, 2534, 4154, 4694,
                                                                       2924, 4739, 1682, 1736,
                                                                       3254, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 5324, 0, 3, 6,
                                                                       4199, 2594, 4289, 4784,
                                                                       2984, 4919, 1844, 1952,
                                                                       3344, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 5594, 0, 3, 6,
                                                                       4289, 2654, 4379, 4919,
                                                                       3074, 5054, 1952, 2060,
                                                                       3524, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 5864, 0, 3, 6,
                                                                       4379, 2714, 4469, 5054,
                                                                       3164, 5189, 2060, 2168,
                                                                       3704, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 6134, 6, 2384,
                                                                       2394, 3914, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 6155, 6, 2394,
                                                                       2404, 3929, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 6176, 6, 2404,
                                                                       2414, 3944, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 6197, 6, 2414,
                                                                       2424, 3959, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 6218, 3, 6, 6134,
                                                                       3914, 6155, 2444, 2474,
                                                                       4064, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 6281, 3, 6, 6155,
                                                                       3929, 6176, 2474, 2504,
                                                                       4109, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 6344, 3, 6, 6176,
                                                                       3944, 6197, 2504, 2534,
                                                                       4154, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 6407, 3, 6, 6218,
                                                                       4064, 6281, 2594, 2654,
                                                                       4379, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 6533, 3, 6, 6281,
                                                                       4109, 6344, 2654, 2714,
                                                                       4469, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 6659, 0, 6, 6134,
                                                                       3914, 6155, 2834, 2864,
                                                                       4649, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 6722, 0, 6, 6155,
                                                                       3929, 6176, 2864, 2894,
                                                                       4694, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 6785, 0, 6, 6176,
                                                                       3944, 6197, 2894, 2924,
                                                                       4739, ncols, gamma, p,
                                                                       q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 6848, 0, 3, 6,
                                                                       6218, 4064, 6281, 6659,
                                                                       4649, 6722, 2984, 3074,
                                                                       5054, ncols, gamma, p,
                                                                       q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 7037, 0, 3, 6,
                                                                       6281, 4109, 6344, 6722,
                                                                       4694, 6785, 3074, 3164,
                                                                       5189, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdh_three_center_electron_repulsion_0(buffer, 7226, 0, 3, 6,
                                                                       6407, 4379, 6533, 6848,
                                                                       5054, 7037, 3344, 3524,
                                                                       5864, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 7604, 6, 3884,
                                                                       3899, 6134, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 7632, 6, 3899,
                                                                       3914, 6155, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 7660, 6, 3914,
                                                                       3929, 6176, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 7688, 6, 3929,
                                                                       3944, 6197, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 7716, 3, 6, 7604,
                                                                       6134, 7632, 3974, 4019,
                                                                       6218, ncols, gamma, p,
                                                                       q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 7800, 3, 6, 7632,
                                                                       6155, 7660, 4019, 4064,
                                                                       6281, ncols, gamma, p,
                                                                       q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 7884, 3, 6, 7660,
                                                                       6176, 7688, 4064, 4109,
                                                                       6344, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 7968, 3, 6, 7716,
                                                                       6218, 7800, 4199, 4289,
                                                                       6407, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 8136, 3, 6, 7800,
                                                                       6281, 7884, 4289, 4379,
                                                                       6533, ncols, gamma, p,
                                                                       q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 8304, 0, 6, 7604,
                                                                       6134, 7632, 4559, 4604,
                                                                       6659, ncols, gamma, p,
                                                                       q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 8388, 0, 6, 7632,
                                                                       6155, 7660, 4604, 4649,
                                                                       6722, ncols, gamma, p,
                                                                       q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 8472, 0, 6, 7660,
                                                                       6176, 7688, 4649, 4694,
                                                                       6785, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 8556, 0, 3, 6,
                                                                       7716, 6218, 7800, 8304,
                                                                       6659, 8388, 4784, 4919,
                                                                       6848, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 8808, 0, 3, 6,
                                                                       7800, 6281, 7884, 8388,
                                                                       6722, 8472, 4919, 5054,
                                                                       7037, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdi_three_center_electron_repulsion_0(buffer, 9060, 0, 3, 6,
                                                                       7968, 6407, 8136, 8556,
                                                                       6848, 8808, 5324, 5594,
                                                                       7226, ncols, gamma, p,
                                                                       q);

                    simdgeo::geom_s_x(buffer, 9564, 9060, 1, 168, ncols, alpha);

                    simdgeo::geom_s_y(buffer, 9732, 9060, 1, 168, ncols, alpha);

                    simdgeo::geom_s_z(buffer, 9900, 9060, 1, 168, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 10068, 9564, 504, ncols);
                }
            }
        }

        simdtrf::transform_i_inner(buffer, 10572, 10068, 6, 1, nmax);

        simdtrf::transform_d_outer(values + n * npairs, nvalues, buffer, 10572, 13, nmax);

        simdtrf::transform_i_inner(buffer, 10572, 10236, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 65 * nvalues + n * npairs, nvalues, buffer, 10572,
                                   13, nmax);

        simdtrf::transform_i_inner(buffer, 10572, 10404, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 130 * nvalues + n * npairs, nvalues, buffer, 10572,
                                   13, nmax);
    }

    for (size_t m = 0; m < 195; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
