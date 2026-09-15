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


#include "SimdThreeCenterElectronRepulsionRecSFH.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

#include "SimdThreeCenterElectronRepulsionVrrRecSDD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformH.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_sfh_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_sfh_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto dimensions = simdfunc::make_column_dimensions(
        a_function, b_function, c_function, npairs, coordinates,
        screenfunc::three_center_electron_repulsion_primitive_bound,
        threshold / static_cast<double>(nprims));

    const auto nmax = simdfunc::prepare_buffer(buffer, 3230, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 77 * natoms * npairs, 0.0);

        return;
    }

    const auto pi = mathconst::pi_value();

    // NOTE: a row of the values spans every atom pair of every atom on the ket
    // side, so a kernel handed the block of one atom steps by this to reach the
    // next component -- which is what lets it be the kernel a two-center form
    // uses, unchanged.

    const auto nvalues = natoms * npairs;

    for (size_t n = 0; n < natoms; n++)
    {
        simdfunc::prepare_buffer(buffer, 3230, 2910, 210, dimensions);

        for (size_t i = 0; i < nprim_a; i++)
        {
            for (size_t j = 0; j < nprim_b; j++)
            {
                const auto p = a_exps[i] + b_exps[j];

                const auto mu = a_exps[i] * b_exps[j] / p;

                const auto fovl = a_norms[i] * b_norms[j];

                const auto fb = a_exps[i] / p;

                const auto fc = b_exps[j] / p;

                simdfunc::compute_pb(buffer, coordinates, 0, nmax, fb);

                simdfunc::compute_pc(buffer, coordinates, c_coordinates, 3, n, nmax, fc);

                simdfunc::compute_pair_exponent(buffer, coordinates, 6, nmax, mu);

                for (size_t k = 0; k < nprim_c; k++)
                {
                    const auto ncols = dimensions[(i * nprim_b + j) * nprim_c + k];

                    if (ncols == 0) continue;

                    const auto gamma = c_exps[k];

                    const auto q = p + gamma;

                    const auto fq = p * gamma / q;

                    const auto fj = 2.0 * fovl * c_norms[k] * pi * pi * std::sqrt(pi)
                                    / (p * gamma * std::sqrt(q));

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 7, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8}, ncols, fj, 6, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 16, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 19, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 22, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 25, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 28, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 31, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 34, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 37, 0, 3, 8, 9,
                                                                       16, 19, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 43, 0, 3, 9, 10,
                                                                       19, 22, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 49, 0, 3, 10, 11,
                                                                       22, 25, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 55, 0, 3, 11, 12,
                                                                       25, 28, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 61, 0, 3, 12, 13,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 67, 0, 3, 13, 14,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 73, 0, 3, 16, 19,
                                                                       37, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 83, 0, 3, 19, 22,
                                                                       43, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 93, 0, 3, 22, 25,
                                                                       49, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 103, 0, 3, 25, 28,
                                                                       55, 61, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 113, 0, 3, 28, 31,
                                                                       61, 67, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 123, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 126, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 129, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 132, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 135, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 138, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 141, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 144, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 147, 3, 10, 22,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 156, 3, 11, 25,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 165, 3, 12, 28,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 174, 3, 13, 31,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 183, 3, 14, 34,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 192, 3, 16, 37,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 210, 3, 19, 43,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 228, 3, 22, 49,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 246, 3, 25, 55,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 264, 3, 28, 61,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 282, 3, 31, 67,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 300, 3, 37, 73,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 330, 3, 43, 83,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 360, 3, 49, 93,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 390, 3, 55, 103,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 420, 3, 61, 113,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 450, 3, 8, 9, 129,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 456, 3, 9, 10,
                                                                       132, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 462, 3, 10, 11,
                                                                       135, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 468, 3, 11, 12,
                                                                       138, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 474, 3, 12, 13,
                                                                       141, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 480, 3, 13, 14,
                                                                       144, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 486, 0, 3, 450,
                                                                       129, 456, 147, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 504, 0, 3, 456,
                                                                       132, 462, 156, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 522, 0, 3, 462,
                                                                       135, 468, 165, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 540, 0, 3, 468,
                                                                       138, 474, 174, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 558, 0, 3, 474,
                                                                       141, 480, 183, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 576, 0, 3, 486,
                                                                       147, 504, 37, 43, 228,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 612, 0, 3, 504,
                                                                       156, 522, 43, 49, 246,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 648, 0, 3, 522,
                                                                       165, 540, 49, 55, 264,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 684, 0, 3, 540,
                                                                       174, 558, 55, 61, 282,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 720, 0, 3, 576,
                                                                       228, 612, 73, 83, 360,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 780, 0, 3, 612,
                                                                       246, 648, 83, 93, 390,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 840, 0, 3, 648,
                                                                       264, 684, 93, 103, 420,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 900, 3, 123, 126,
                                                                       450, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 910, 3, 126, 129,
                                                                       456, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 920, 3, 129, 132,
                                                                       462, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 930, 3, 132, 135,
                                                                       468, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 940, 3, 135, 138,
                                                                       474, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 950, 3, 138, 141,
                                                                       480, ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 960, 0, 3, 900,
                                                                       450, 910, 486, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 990, 0, 3, 910,
                                                                       456, 920, 504, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 1020, 0, 3, 920,
                                                                       462, 930, 522, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 1050, 0, 3, 930,
                                                                       468, 940, 540, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 1080, 0, 3, 940,
                                                                       474, 950, 558, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 1110, 0, 3, 960,
                                                                       486, 990, 192, 210, 576,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 1170, 0, 3, 990,
                                                                       504, 1020, 210, 228, 612,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 1230, 0, 3, 1020,
                                                                       522, 1050, 228, 246, 648,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 1290, 0, 3, 1050,
                                                                       540, 1080, 246, 264, 684,
                                                                       ncols, gamma, p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 1350, 0, 3, 1110,
                                                                       576, 1170, 300, 330, 720,
                                                                       ncols, gamma, p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 1450, 0, 3, 1170,
                                                                       612, 1230, 330, 360, 780,
                                                                       ncols, gamma, p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 1550, 0, 3, 1230,
                                                                       648, 1290, 360, 390, 840,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 1650, 3, 450, 456,
                                                                       920, ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 1665, 3, 456, 462,
                                                                       930, ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 1680, 3, 462, 468,
                                                                       940, ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 1695, 3, 468, 474,
                                                                       950, ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 1710, 0, 3, 1650,
                                                                       920, 1665, 486, 504, 1020,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 1755, 0, 3, 1665,
                                                                       930, 1680, 504, 522, 1050,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 1800, 0, 3, 1680,
                                                                       940, 1695, 522, 540, 1080,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 1845, 0, 3, 1710,
                                                                       1020, 1755, 576, 612,
                                                                       1230, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 1935, 0, 3, 1755,
                                                                       1050, 1800, 612, 648,
                                                                       1290, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 2025, 0, 3, 1845,
                                                                       1230, 1935, 720, 780,
                                                                       1550, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 2175, 3, 900, 910,
                                                                       1650, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 2196, 3, 910, 920,
                                                                       1665, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 2217, 3, 920, 930,
                                                                       1680, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 2238, 3, 930, 940,
                                                                       1695, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 2259, 0, 3, 2175,
                                                                       1650, 2196, 960, 990,
                                                                       1710, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 2322, 0, 3, 2196,
                                                                       1665, 2217, 990, 1020,
                                                                       1755, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 2385, 0, 3, 2217,
                                                                       1680, 2238, 1020, 1050,
                                                                       1800, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 2448, 0, 3, 2259,
                                                                       1710, 2322, 1110, 1170,
                                                                       1845, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 2574, 0, 3, 2322,
                                                                       1755, 2385, 1170, 1230,
                                                                       1935, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 2700, 0, 3, 2448,
                                                                       1845, 2574, 1350, 1450,
                                                                       2025, ncols, gamma, p,
                                                                       q);

                    simdfunc::contract_primitives(buffer, 2910, 2700, 210, ncols);
                }
            }
        }

        simdtrf::transform_h_inner(buffer, 3120, 2910, 10, 1, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 3120, 11, nmax);
    }

    for (size_t m = 0; m < 77; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
