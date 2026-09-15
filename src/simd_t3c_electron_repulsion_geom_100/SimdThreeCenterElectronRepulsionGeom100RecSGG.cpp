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


#include "SimdThreeCenterElectronRepulsionGeom100RecSGG.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecPFD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPFF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPFG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPGD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPGF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPGG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPGP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPGS.hpp"
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
#include "SimdThreeCenterElectronRepulsionVrrRecSFD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransformG.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_geom_100_sgg_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_geom_100_sgg_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 15594, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 243 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 15594, 14784, 675, dimensions);

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

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 95, 3, 6, 20, 23,
                                                                       47, 53, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 105, 3, 6, 23, 26,
                                                                       53, 59, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 115, 3, 6, 26, 29,
                                                                       59, 65, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 125, 3, 6, 29, 32,
                                                                       65, 71, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 135, 3, 6, 32, 35,
                                                                       71, 77, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 145, 3, 6, 35, 38,
                                                                       77, 83, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 155, 3, 6, 38, 41,
                                                                       83, 89, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 165, 3, 6, 47, 53,
                                                                       95, 105, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 180, 3, 6, 53, 59,
                                                                       105, 115, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 195, 3, 6, 59, 65,
                                                                       115, 125, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 210, 3, 6, 65, 71,
                                                                       125, 135, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 225, 3, 6, 71, 77,
                                                                       135, 145, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 240, 3, 6, 77, 83,
                                                                       145, 155, ncols, gamma, p,
                                                                       q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 255, 0, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 258, 0, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 261, 0, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 264, 0, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 267, 0, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 270, 0, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 273, 0, 6, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 276, 0, 6, 10, 11,
                                                                       20, 23, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 285, 0, 6, 11, 12,
                                                                       23, 26, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 294, 0, 6, 12, 13,
                                                                       26, 29, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 303, 0, 6, 13, 14,
                                                                       29, 32, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 312, 0, 6, 14, 15,
                                                                       32, 35, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 321, 0, 6, 15, 16,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 330, 0, 6, 16, 17,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 339, 0, 6, 17, 18,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 348, 0, 3, 6, 20,
                                                                       23, 47, 53, 276, 285,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 366, 0, 3, 6, 23,
                                                                       26, 53, 59, 285, 294,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 384, 0, 3, 6, 26,
                                                                       29, 59, 65, 294, 303,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 402, 0, 3, 6, 29,
                                                                       32, 65, 71, 303, 312,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 420, 0, 3, 6, 32,
                                                                       35, 71, 77, 312, 321,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 438, 0, 3, 6, 35,
                                                                       38, 77, 83, 321, 330,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 456, 0, 3, 6, 38,
                                                                       41, 83, 89, 330, 339,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 474, 0, 3, 6, 47,
                                                                       53, 95, 105, 348, 366,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 504, 0, 3, 6, 53,
                                                                       59, 105, 115, 366, 384,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 534, 0, 3, 6, 59,
                                                                       65, 115, 125, 384, 402,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 564, 0, 3, 6, 65,
                                                                       71, 125, 135, 402, 420,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 594, 0, 3, 6, 71,
                                                                       77, 135, 145, 420, 438,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 624, 0, 3, 6, 77,
                                                                       83, 145, 155, 438, 456,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 654, 0, 3, 6, 95,
                                                                       105, 165, 180, 474, 504,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 699, 0, 3, 6, 105,
                                                                       115, 180, 195, 504, 534,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 744, 0, 3, 6, 115,
                                                                       125, 195, 210, 534, 564,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 789, 0, 3, 6, 125,
                                                                       135, 210, 225, 564, 594,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 834, 0, 3, 6, 135,
                                                                       145, 225, 240, 594, 624,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 879, 6, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 882, 6, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 885, 6, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 888, 6, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 891, 6, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 894, 6, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 897, 6, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 900, 6, 19, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 903, 6, 12, 26,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 912, 6, 13, 29,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 921, 6, 14, 32,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 930, 6, 15, 35,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 939, 6, 16, 38,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 948, 6, 17, 41,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 957, 6, 18, 44,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 966, 6, 26, 59,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 984, 6, 29, 65,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1002, 6, 32, 71,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1020, 6, 35, 77,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1038, 6, 38, 83,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1056, 6, 41, 89,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1074, 6, 59, 115,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1104, 6, 65, 125,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1134, 6, 71, 135,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1164, 6, 77, 145,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1194, 6, 83, 155,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1224, 6, 115, 195,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1269, 6, 125, 210,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1314, 6, 135, 225,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1359, 6, 145, 240,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1404, 6, 12, 255,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1413, 6, 13, 258,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1422, 6, 14, 261,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1431, 6, 15, 264,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1440, 6, 16, 267,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1449, 6, 17, 270,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1458, 6, 18, 273,
                                                                       ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1467, 6, 26, 255,
                                                                       294, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1494, 6, 29, 258,
                                                                       303, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1521, 6, 32, 261,
                                                                       312, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1548, 6, 35, 264,
                                                                       321, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1575, 6, 38, 267,
                                                                       330, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1602, 6, 41, 270,
                                                                       339, ncols, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1629, 3, 6, 59,
                                                                       1467, 294, 1494, 384,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1683, 3, 6, 65,
                                                                       1494, 303, 1521, 402,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1737, 3, 6, 71,
                                                                       1521, 312, 1548, 420,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1791, 3, 6, 77,
                                                                       1548, 321, 1575, 438,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1845, 3, 6, 83,
                                                                       1575, 330, 1602, 456,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 1899, 3, 6, 115,
                                                                       1629, 384, 1683, 534,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 1989, 3, 6, 125,
                                                                       1683, 402, 1737, 564,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 2079, 3, 6, 135,
                                                                       1737, 420, 1791, 594,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 2169, 3, 6, 145,
                                                                       1791, 438, 1845, 624,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 2259, 3, 6, 195,
                                                                       1899, 534, 1989, 744,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 2394, 3, 6, 210,
                                                                       1989, 564, 2079, 789,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 2529, 3, 6, 225,
                                                                       2079, 594, 2169, 834,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2664, 6, 10, 11,
                                                                       879, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2670, 6, 11, 12,
                                                                       882, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2676, 6, 12, 13,
                                                                       885, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2682, 6, 13, 14,
                                                                       888, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2688, 6, 14, 15,
                                                                       891, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2694, 6, 15, 16,
                                                                       894, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2700, 6, 16, 17,
                                                                       897, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2706, 6, 17, 18,
                                                                       900, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2712, 3, 6, 2664,
                                                                       879, 2670, 903, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2730, 3, 6, 2670,
                                                                       882, 2676, 912, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2748, 3, 6, 2676,
                                                                       885, 2682, 921, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2766, 3, 6, 2682,
                                                                       888, 2688, 930, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2784, 3, 6, 2688,
                                                                       891, 2694, 939, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2802, 3, 6, 2694,
                                                                       894, 2700, 948, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2820, 3, 6, 2700,
                                                                       897, 2706, 957, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2838, 3, 6, 2712,
                                                                       903, 2730, 47, 53, 966,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2874, 3, 6, 2730,
                                                                       912, 2748, 53, 59, 984,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2910, 3, 6, 2748,
                                                                       921, 2766, 59, 65, 1002,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2946, 3, 6, 2766,
                                                                       930, 2784, 65, 71, 1020,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2982, 3, 6, 2784,
                                                                       939, 2802, 71, 77, 1038,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3018, 3, 6, 2802,
                                                                       948, 2820, 77, 83, 1056,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3054, 3, 6, 2838,
                                                                       966, 2874, 95, 105, 1074,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3114, 3, 6, 2874,
                                                                       984, 2910, 105, 115, 1104,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3174, 3, 6, 2910,
                                                                       1002, 2946, 115, 125,
                                                                       1134, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3234, 3, 6, 2946,
                                                                       1020, 2982, 125, 135,
                                                                       1164, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3294, 3, 6, 2982,
                                                                       1038, 3018, 135, 145,
                                                                       1194, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3354, 3, 6, 3054,
                                                                       1074, 3114, 165, 180,
                                                                       1224, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3444, 3, 6, 3114,
                                                                       1104, 3174, 180, 195,
                                                                       1269, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3534, 3, 6, 3174,
                                                                       1134, 3234, 195, 210,
                                                                       1314, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3624, 3, 6, 3234,
                                                                       1164, 3294, 210, 225,
                                                                       1359, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3714, 0, 6, 2664,
                                                                       879, 2670, 1404, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3732, 0, 6, 2670,
                                                                       882, 2676, 1413, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3750, 0, 6, 2676,
                                                                       885, 2682, 1422, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3768, 0, 6, 2682,
                                                                       888, 2688, 1431, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3786, 0, 6, 2688,
                                                                       891, 2694, 1440, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3804, 0, 6, 2694,
                                                                       894, 2700, 1449, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3822, 0, 6, 2700,
                                                                       897, 2706, 1458, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3840, 0, 3, 6,
                                                                       2712, 903, 2730, 3714,
                                                                       1404, 3732, 276, 285,
                                                                       1467, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3894, 0, 3, 6,
                                                                       2730, 912, 2748, 3732,
                                                                       1413, 3750, 285, 294,
                                                                       1494, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3948, 0, 3, 6,
                                                                       2748, 921, 2766, 3750,
                                                                       1422, 3768, 294, 303,
                                                                       1521, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 4002, 0, 3, 6,
                                                                       2766, 930, 2784, 3768,
                                                                       1431, 3786, 303, 312,
                                                                       1548, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 4056, 0, 3, 6,
                                                                       2784, 939, 2802, 3786,
                                                                       1440, 3804, 312, 321,
                                                                       1575, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 4110, 0, 3, 6,
                                                                       2802, 948, 2820, 3804,
                                                                       1449, 3822, 321, 330,
                                                                       1602, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 4164, 0, 3, 6,
                                                                       2838, 966, 2874, 3840,
                                                                       1467, 3894, 348, 366,
                                                                       1629, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 4272, 0, 3, 6,
                                                                       2874, 984, 2910, 3894,
                                                                       1494, 3948, 366, 384,
                                                                       1683, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 4380, 0, 3, 6,
                                                                       2910, 1002, 2946, 3948,
                                                                       1521, 4002, 384, 402,
                                                                       1737, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 4488, 0, 3, 6,
                                                                       2946, 1020, 2982, 4002,
                                                                       1548, 4056, 402, 420,
                                                                       1791, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 4596, 0, 3, 6,
                                                                       2982, 1038, 3018, 4056,
                                                                       1575, 4110, 420, 438,
                                                                       1845, ncols, gamma, p,
                                                                       q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 4704, 0, 3, 6,
                                                                       3054, 1074, 3114, 3840,
                                                                       3894, 4164, 1629, 4272,
                                                                       474, 504, 1899, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 4884, 0, 3, 6,
                                                                       3114, 1104, 3174, 3894,
                                                                       3948, 4272, 1683, 4380,
                                                                       504, 534, 1989, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 5064, 0, 3, 6,
                                                                       3174, 1134, 3234, 3948,
                                                                       4002, 4380, 1737, 4488,
                                                                       534, 564, 2079, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 5244, 0, 3, 6,
                                                                       3234, 1164, 3294, 4002,
                                                                       4056, 4488, 1791, 4596,
                                                                       564, 594, 2169, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgd_three_center_electron_repulsion_0(buffer, 5424, 0, 3, 6,
                                                                       3354, 1224, 3444, 4164,
                                                                       4272, 4704, 1899, 4884,
                                                                       654, 699, 2259, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgd_three_center_electron_repulsion_0(buffer, 5694, 0, 3, 6,
                                                                       3444, 1269, 3534, 4272,
                                                                       4380, 4884, 1989, 5064,
                                                                       699, 744, 2394, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgd_three_center_electron_repulsion_0(buffer, 5964, 0, 3, 6,
                                                                       3534, 1314, 3624, 4380,
                                                                       4488, 5064, 2079, 5244,
                                                                       744, 789, 2529, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6234, 6, 879, 882,
                                                                       2676, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6244, 6, 882, 885,
                                                                       2682, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6254, 6, 885, 888,
                                                                       2688, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6264, 6, 888, 891,
                                                                       2694, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6274, 6, 891, 894,
                                                                       2700, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6284, 6, 894, 897,
                                                                       2706, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6294, 3, 6, 6234,
                                                                       2676, 6244, 2748, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6324, 3, 6, 6244,
                                                                       2682, 6254, 2766, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6354, 3, 6, 6254,
                                                                       2688, 6264, 2784, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6384, 3, 6, 6264,
                                                                       2694, 6274, 2802, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6414, 3, 6, 6274,
                                                                       2700, 6284, 2820, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 6444, 3, 6, 6294,
                                                                       2748, 6324, 966, 984,
                                                                       2910, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 6504, 3, 6, 6324,
                                                                       2766, 6354, 984, 1002,
                                                                       2946, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 6564, 3, 6, 6354,
                                                                       2784, 6384, 1002, 1020,
                                                                       2982, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 6624, 3, 6, 6384,
                                                                       2802, 6414, 1020, 1038,
                                                                       3018, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 6684, 3, 6, 6444,
                                                                       2910, 6504, 1074, 1104,
                                                                       3174, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 6784, 3, 6, 6504,
                                                                       2946, 6564, 1104, 1134,
                                                                       3234, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 6884, 3, 6, 6564,
                                                                       2982, 6624, 1134, 1164,
                                                                       3294, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 6984, 3, 6, 6684,
                                                                       3174, 6784, 1224, 1269,
                                                                       3534, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 7134, 3, 6, 6784,
                                                                       3234, 6884, 1269, 1314,
                                                                       3624, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 7284, 0, 6, 6234,
                                                                       2676, 6244, 1404, 1413,
                                                                       3750, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 7314, 0, 6, 6244,
                                                                       2682, 6254, 1413, 1422,
                                                                       3768, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 7344, 0, 6, 6254,
                                                                       2688, 6264, 1422, 1431,
                                                                       3786, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 7374, 0, 6, 6264,
                                                                       2694, 6274, 1431, 1440,
                                                                       3804, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 7404, 0, 6, 6274,
                                                                       2700, 6284, 1440, 1449,
                                                                       3822, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 7434, 0, 3, 6,
                                                                       6294, 2748, 6324, 7284,
                                                                       3750, 7314, 1467, 1494,
                                                                       3948, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 7524, 0, 3, 6,
                                                                       6324, 2766, 6354, 7314,
                                                                       3768, 7344, 1494, 1521,
                                                                       4002, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 7614, 0, 3, 6,
                                                                       6354, 2784, 6384, 7344,
                                                                       3786, 7374, 1521, 1548,
                                                                       4056, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 7704, 0, 3, 6,
                                                                       6384, 2802, 6414, 7374,
                                                                       3804, 7404, 1548, 1575,
                                                                       4110, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 7794, 0, 3, 6,
                                                                       6444, 2910, 6504, 7434,
                                                                       3948, 7524, 1629, 1683,
                                                                       4380, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 7974, 0, 3, 6,
                                                                       6504, 2946, 6564, 7524,
                                                                       4002, 7614, 1683, 1737,
                                                                       4488, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 8154, 0, 3, 6,
                                                                       6564, 2982, 6624, 7614,
                                                                       4056, 7704, 1737, 1791,
                                                                       4596, ncols, gamma, p,
                                                                       q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 8334, 0, 3, 6,
                                                                       6684, 3174, 6784, 7434,
                                                                       7524, 7794, 4380, 7974,
                                                                       1899, 1989, 5064, ncols,
                                                                       gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 8634, 0, 3, 6,
                                                                       6784, 3234, 6884, 7524,
                                                                       7614, 7974, 4488, 8154,
                                                                       1989, 2079, 5244, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgf_three_center_electron_repulsion_0(buffer, 8934, 0, 3, 6,
                                                                       6984, 3534, 7134, 7794,
                                                                       7974, 8334, 5064, 8634,
                                                                       2259, 2394, 5964, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9384, 6, 2664,
                                                                       2670, 6234, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9399, 6, 2670,
                                                                       2676, 6244, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9414, 6, 2676,
                                                                       2682, 6254, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9429, 6, 2682,
                                                                       2688, 6264, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9444, 6, 2688,
                                                                       2694, 6274, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9459, 6, 2694,
                                                                       2700, 6284, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 9474, 3, 6, 9384,
                                                                       6234, 9399, 2712, 2730,
                                                                       6294, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 9519, 3, 6, 9399,
                                                                       6244, 9414, 2730, 2748,
                                                                       6324, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 9564, 3, 6, 9414,
                                                                       6254, 9429, 2748, 2766,
                                                                       6354, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 9609, 3, 6, 9429,
                                                                       6264, 9444, 2766, 2784,
                                                                       6384, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 9654, 3, 6, 9444,
                                                                       6274, 9459, 2784, 2802,
                                                                       6414, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 9699, 3, 6, 9474,
                                                                       6294, 9519, 2838, 2874,
                                                                       6444, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 9789, 3, 6, 9519,
                                                                       6324, 9564, 2874, 2910,
                                                                       6504, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 9879, 3, 6, 9564,
                                                                       6354, 9609, 2910, 2946,
                                                                       6564, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 9969, 3, 6, 9609,
                                                                       6384, 9654, 2946, 2982,
                                                                       6624, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 10059, 3, 6, 9699,
                                                                       6444, 9789, 3054, 3114,
                                                                       6684, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 10209, 3, 6, 9789,
                                                                       6504, 9879, 3114, 3174,
                                                                       6784, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 10359, 3, 6, 9879,
                                                                       6564, 9969, 3174, 3234,
                                                                       6884, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 10509, 3, 6,
                                                                       10059, 6684, 10209, 3354,
                                                                       3444, 6984, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 10734, 3, 6,
                                                                       10209, 6784, 10359, 3444,
                                                                       3534, 7134, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 10959, 0, 6, 9384,
                                                                       6234, 9399, 3714, 3732,
                                                                       7284, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 11004, 0, 6, 9399,
                                                                       6244, 9414, 3732, 3750,
                                                                       7314, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 11049, 0, 6, 9414,
                                                                       6254, 9429, 3750, 3768,
                                                                       7344, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 11094, 0, 6, 9429,
                                                                       6264, 9444, 3768, 3786,
                                                                       7374, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 11139, 0, 6, 9444,
                                                                       6274, 9459, 3786, 3804,
                                                                       7404, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 11184, 0, 3, 6,
                                                                       9474, 6294, 9519, 10959,
                                                                       7284, 11004, 3840, 3894,
                                                                       7434, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 11319, 0, 3, 6,
                                                                       9519, 6324, 9564, 11004,
                                                                       7314, 11049, 3894, 3948,
                                                                       7524, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 11454, 0, 3, 6,
                                                                       9564, 6354, 9609, 11049,
                                                                       7344, 11094, 3948, 4002,
                                                                       7614, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 11589, 0, 3, 6,
                                                                       9609, 6384, 9654, 11094,
                                                                       7374, 11139, 4002, 4056,
                                                                       7704, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 11724, 0, 3, 6,
                                                                       9699, 6444, 9789, 11184,
                                                                       7434, 11319, 4164, 4272,
                                                                       7794, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 11994, 0, 3, 6,
                                                                       9789, 6504, 9879, 11319,
                                                                       7524, 11454, 4272, 4380,
                                                                       7974, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 12264, 0, 3, 6,
                                                                       9879, 6564, 9969, 11454,
                                                                       7614, 11589, 4380, 4488,
                                                                       8154, ncols, gamma, p,
                                                                       q);

                    compute_prim_pfg_three_center_electron_repulsion_0(buffer, 12534, 0, 3, 6,
                                                                       10059, 6684, 10209, 11184,
                                                                       11319, 11724, 7794, 11994,
                                                                       4704, 4884, 8334, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfg_three_center_electron_repulsion_0(buffer, 12984, 0, 3, 6,
                                                                       10209, 6784, 10359, 11319,
                                                                       11454, 11994, 7974, 12264,
                                                                       4884, 5064, 8634, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgg_three_center_electron_repulsion_0(buffer, 13434, 0, 3, 6,
                                                                       10509, 6984, 10734, 11724,
                                                                       11994, 12534, 8334, 12984,
                                                                       5424, 5694, 8934, ncols,
                                                                       gamma, p, q);

                    simdgeo::geom_s_x(buffer, 14109, 13434, 1, 225, ncols, alpha);

                    simdgeo::geom_s_y(buffer, 14334, 13434, 1, 225, ncols, alpha);

                    simdgeo::geom_s_z(buffer, 14559, 13434, 1, 225, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 14784, 14109, 675, ncols);
                }
            }
        }

        simdtrf::transform_g_inner(buffer, 15459, 14784, 15, 1, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 15459, 9, nmax);

        simdtrf::transform_g_inner(buffer, 15459, 15009, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 81 * nvalues + n * npairs, nvalues, buffer, 15459, 9,
                                   nmax);

        simdtrf::transform_g_inner(buffer, 15459, 15234, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 162 * nvalues + n * npairs, nvalues, buffer, 15459,
                                   9, nmax);
    }

    for (size_t m = 0; m < 243; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
