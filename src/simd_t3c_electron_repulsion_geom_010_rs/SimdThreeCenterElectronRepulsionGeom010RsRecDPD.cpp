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


#include "SimdThreeCenterElectronRepulsionGeom010RsRecDPD.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

#include "SimdGeometryD1.hpp"
#include "SimdGeometryF1.hpp"
#include "SimdGeometryP1.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferGeom010XDP.hpp"
#include "SimdTransferGeom010XPD.hpp"
#include "SimdTransferGeom010XPP.hpp"
#include "SimdTransferGeom010YDP.hpp"
#include "SimdTransferGeom010YPD.hpp"
#include "SimdTransferGeom010YPP.hpp"
#include "SimdTransferGeom010ZDP.hpp"
#include "SimdTransferGeom010ZPD.hpp"
#include "SimdTransferGeom010ZPP.hpp"
#include "SimdTransferPP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_010_dpd_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_010_dpd_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 5214, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 450 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 5214, 2232, 1402, dimensions);

        for (size_t i = 0; i < nprim_a; i++)
        {
            for (size_t j = 0; j < nprim_b; j++)
            {
                const auto p = a_exps[i] + b_exps[j];

                const auto beta = b_exps[j];

                const auto fovl = a_norms[i] * b_norms[j];

                const auto fb = a_exps[i] / p;

                const auto fc = b_exps[j] / p;

                simdfunc::compute_pb(buffer, coordinates, 0, nmax, fb);

                simdfunc::compute_pc(buffer, coordinates, c_coordinates, 3, n, nmax, fc);

                for (size_t k = 0; k < nprim_c; k++)
                {
                    const auto ncols = dimensions[(i * nprim_b + j) * nprim_c + k];

                    if (ncols == 0) continue;

                    const auto gamma = c_exps[k];

                    const auto q = p + gamma;

                    const auto fq = p * gamma / q;

                    const auto fj = 2.0 * fovl * c_norms[k] * pi * pi * std::sqrt(pi)
                                    / (p * gamma * std::sqrt(q));

                    simdfunc::compute_full_t3c_erf_boys_function(buffer, coordinates, 6, 3, 6,
                                                                 ncols, fj, i * nprim_b + j, fq,
                                                                 omega);

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 14, 3, 6,
                                                             ncols, fj, i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 22, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 25, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 28, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 31, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 34, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 37, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 40, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 43, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 46, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 49, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 52, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 55, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 58, 0, 3, 7, 8,
                                                                       22, 25, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 64, 0, 3, 8, 9,
                                                                       25, 28, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 70, 0, 3, 9, 10,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 76, 0, 3, 10, 11,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 82, 0, 3, 11, 12,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 88, 0, 3, 15, 16,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 94, 0, 3, 16, 17,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 100, 0, 3, 17, 18,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 106, 0, 3, 18, 19,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 112, 0, 3, 19, 20,
                                                                       52, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 118, 0, 3, 22, 25,
                                                                       58, 64, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 128, 0, 3, 25, 28,
                                                                       64, 70, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 138, 0, 3, 28, 31,
                                                                       70, 76, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 148, 0, 3, 31, 34,
                                                                       76, 82, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 158, 0, 3, 40, 43,
                                                                       88, 94, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 168, 0, 3, 43, 46,
                                                                       94, 100, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 178, 0, 3, 46, 49,
                                                                       100, 106, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 188, 0, 3, 49, 52,
                                                                       106, 112, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 198, 0, 3, 58, 64,
                                                                       118, 128, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 213, 0, 3, 64, 70,
                                                                       128, 138, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 228, 0, 3, 70, 76,
                                                                       138, 148, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 243, 0, 3, 88, 94,
                                                                       158, 168, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 258, 0, 3, 94,
                                                                       100, 168, 178, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 273, 0, 3, 100,
                                                                       106, 178, 188, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 288, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 291, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 294, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 297, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 300, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 303, 3, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 306, 3, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 309, 3, 19, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 312, 3, 20, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 315, 3, 21, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 318, 3, 9, 28,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 327, 3, 10, 31,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 336, 3, 11, 34,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 345, 3, 12, 37,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 354, 3, 17, 46,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 363, 3, 18, 49,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 372, 3, 19, 52,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 381, 3, 20, 55,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 390, 3, 28, 70,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 408, 3, 31, 76,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 426, 3, 34, 82,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 444, 3, 46, 100,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 462, 3, 49, 106,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 480, 3, 52, 112,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 498, 3, 70, 138,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 528, 3, 76, 148,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 558, 3, 100, 178,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 588, 3, 106, 188,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 618, 3, 138, 228,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 663, 3, 178, 273,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 708, 3, 7, 8, 288,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 714, 3, 8, 9, 291,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 720, 3, 9, 10,
                                                                       294, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 726, 3, 10, 11,
                                                                       297, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 732, 3, 11, 12,
                                                                       300, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 738, 3, 15, 16,
                                                                       303, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 744, 3, 16, 17,
                                                                       306, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 750, 3, 17, 18,
                                                                       309, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 756, 3, 18, 19,
                                                                       312, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 762, 3, 19, 20,
                                                                       315, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 768, 0, 3, 708,
                                                                       288, 714, 318, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 786, 0, 3, 714,
                                                                       291, 720, 327, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 804, 0, 3, 720,
                                                                       294, 726, 336, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 822, 0, 3, 726,
                                                                       297, 732, 345, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 840, 0, 3, 738,
                                                                       303, 744, 354, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 858, 0, 3, 744,
                                                                       306, 750, 363, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 876, 0, 3, 750,
                                                                       309, 756, 372, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 894, 0, 3, 756,
                                                                       312, 762, 381, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 912, 0, 3, 768,
                                                                       318, 786, 58, 64, 390,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 948, 0, 3, 786,
                                                                       327, 804, 64, 70, 408,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 984, 0, 3, 804,
                                                                       336, 822, 70, 76, 426,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1020, 0, 3, 840,
                                                                       354, 858, 88, 94, 444,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1056, 0, 3, 858,
                                                                       363, 876, 94, 100, 462,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1092, 0, 3, 876,
                                                                       372, 894, 100, 106, 480,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 1128, 0, 3, 912,
                                                                       390, 948, 118, 128, 498,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 1188, 0, 3, 948,
                                                                       408, 984, 128, 138, 528,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 1248, 0, 3, 1020,
                                                                       444, 1056, 158, 168, 558,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 1308, 0, 3, 1056,
                                                                       462, 1092, 168, 178, 588,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 1368, 0, 3, 1128,
                                                                       498, 1188, 198, 213, 618,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 1458, 0, 3, 1248,
                                                                       558, 1308, 243, 258, 663,
                                                                       ncols, gamma, p, q);

                    simdgeo::geom_p_x(buffer, 1548, 708, 912, 1, 6, ncols, beta);

                    simdgeo::geom_p_y(buffer, 1566, 708, 912, 1, 6, ncols, beta);

                    simdgeo::geom_p_z(buffer, 1584, 708, 912, 1, 6, ncols, beta);

                    simdgeo::geom_p_x(buffer, 1602, 738, 1020, 1, 6, ncols, beta);

                    simdgeo::geom_p_y(buffer, 1620, 738, 1020, 1, 6, ncols, beta);

                    simdgeo::geom_p_z(buffer, 1638, 738, 1020, 1, 6, ncols, beta);

                    simdgeo::geom_d_x(buffer, 1656, 768, 1128, 1, 6, ncols, beta);

                    simdgeo::geom_d_y(buffer, 1692, 768, 1128, 1, 6, ncols, beta);

                    simdgeo::geom_d_z(buffer, 1728, 768, 1128, 1, 6, ncols, beta);

                    simdgeo::geom_d_x(buffer, 1764, 840, 1248, 1, 6, ncols, beta);

                    simdgeo::geom_d_y(buffer, 1800, 840, 1248, 1, 6, ncols, beta);

                    simdgeo::geom_d_z(buffer, 1836, 840, 1248, 1, 6, ncols, beta);

                    simdgeo::geom_f_x(buffer, 1872, 912, 1368, 1, 6, ncols, beta);

                    simdgeo::geom_f_y(buffer, 1932, 912, 1368, 1, 6, ncols, beta);

                    simdgeo::geom_f_z(buffer, 1992, 912, 1368, 1, 6, ncols, beta);

                    simdgeo::geom_f_x(buffer, 2052, 1020, 1458, 1, 6, ncols, beta);

                    simdgeo::geom_f_y(buffer, 2112, 1020, 1458, 1, 6, ncols, beta);

                    simdgeo::geom_f_z(buffer, 2172, 1020, 1458, 1, 6, ncols, beta);

                    simdfunc::contract_primitives(buffer, 2232, 1548, 18, ncols);

                    simdfunc::contract_primitives(buffer, 2265, 1566, 18, ncols);

                    simdfunc::contract_primitives(buffer, 2298, 1584, 18, ncols);

                    simdfunc::contract_primitives(buffer, 2331, 768, 18, ncols);

                    simdfunc::contract_primitives(buffer, 2364, 1602, 18, ncols);

                    simdfunc::contract_primitives(buffer, 2397, 1620, 18, ncols);

                    simdfunc::contract_primitives(buffer, 2430, 1638, 18, ncols);

                    simdfunc::contract_primitives(buffer, 2463, 840, 18, ncols);

                    simdfunc::contract_primitives(buffer, 2496, 1656, 36, ncols);

                    simdfunc::contract_primitives(buffer, 2562, 1692, 36, ncols);

                    simdfunc::contract_primitives(buffer, 2628, 1728, 36, ncols);

                    simdfunc::contract_primitives(buffer, 2694, 912, 36, ncols);

                    simdfunc::contract_primitives(buffer, 2760, 1764, 36, ncols);

                    simdfunc::contract_primitives(buffer, 2826, 1800, 36, ncols);

                    simdfunc::contract_primitives(buffer, 2892, 1836, 36, ncols);

                    simdfunc::contract_primitives(buffer, 2958, 1020, 36, ncols);

                    simdfunc::contract_primitives(buffer, 3024, 1872, 60, ncols);

                    simdfunc::contract_primitives(buffer, 3134, 1932, 60, ncols);

                    simdfunc::contract_primitives(buffer, 3244, 1992, 60, ncols);

                    simdfunc::contract_primitives(buffer, 3354, 2052, 60, ncols);

                    simdfunc::contract_primitives(buffer, 3464, 2112, 60, ncols);

                    simdfunc::contract_primitives(buffer, 3574, 2172, 60, ncols);
                }
            }
        }

        simdtrf::transform_d_inner(buffer, 2250, 2232, 3, 1, nmax);

        simdtrf::transform_d_inner(buffer, 2283, 2265, 3, 1, nmax);

        simdtrf::transform_d_inner(buffer, 2316, 2298, 3, 1, nmax);

        simdtrf::transform_d_inner(buffer, 2349, 2331, 3, 1, nmax);

        simdtrf::transform_d_inner(buffer, 2382, 2364, 3, 1, nmax);

        simdtrf::transform_d_inner(buffer, 2415, 2397, 3, 1, nmax);

        simdtrf::transform_d_inner(buffer, 2448, 2430, 3, 1, nmax);

        simdtrf::transform_d_inner(buffer, 2481, 2463, 3, 1, nmax);

        simdtrf::transform_d_inner(buffer, 2532, 2496, 6, 1, nmax);

        simdtrf::transform_d_inner(buffer, 2598, 2562, 6, 1, nmax);

        simdtrf::transform_d_inner(buffer, 2664, 2628, 6, 1, nmax);

        simdtrf::transform_d_inner(buffer, 2730, 2694, 6, 1, nmax);

        simdtrf::transform_d_inner(buffer, 2796, 2760, 6, 1, nmax);

        simdtrf::transform_d_inner(buffer, 2862, 2826, 6, 1, nmax);

        simdtrf::transform_d_inner(buffer, 2928, 2892, 6, 1, nmax);

        simdtrf::transform_d_inner(buffer, 2994, 2958, 6, 1, nmax);

        simdtrf::transform_d_inner(buffer, 3084, 3024, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 3194, 3134, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 3304, 3244, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 3414, 3354, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 3524, 3464, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 3634, 3574, 10, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pp(buffer, coordinates, 3684, 2250, 2349, 2532, 5, nmax);

        simdtrf::compute_hrr_geom_010y_pp(buffer, coordinates, 3729, 2283, 2349, 2598, 5, nmax);

        simdtrf::compute_hrr_geom_010z_pp(buffer, coordinates, 3774, 2316, 2349, 2664, 5, nmax);

        simdtrf::compute_hrr_pp(buffer, coordinates, 3819, 2349, 2730, 5, nmax);

        simdtrf::compute_hrr_geom_010x_pp(buffer, coordinates, 3864, 2382, 2481, 2796, 5, nmax);

        simdtrf::compute_hrr_geom_010y_pp(buffer, coordinates, 3909, 2415, 2481, 2862, 5, nmax);

        simdtrf::compute_hrr_geom_010z_pp(buffer, coordinates, 3954, 2448, 2481, 2928, 5, nmax);

        simdtrf::compute_hrr_pp(buffer, coordinates, 3999, 2481, 2994, 5, nmax);

        simdtrf::compute_hrr_geom_010x_pd(buffer, coordinates, 4044, 2532, 2730, 3084, 5, nmax);

        simdtrf::compute_hrr_geom_010y_pd(buffer, coordinates, 4134, 2598, 2730, 3194, 5, nmax);

        simdtrf::compute_hrr_geom_010z_pd(buffer, coordinates, 4224, 2664, 2730, 3304, 5, nmax);

        simdtrf::compute_hrr_geom_010x_pd(buffer, coordinates, 4314, 2796, 2994, 3414, 5, nmax);

        simdtrf::compute_hrr_geom_010y_pd(buffer, coordinates, 4404, 2862, 2994, 3524, 5, nmax);

        simdtrf::compute_hrr_geom_010z_pd(buffer, coordinates, 4494, 2928, 2994, 3634, 5, nmax);

        simdtrf::compute_hrr_geom_010x_dp_out_of_second(buffer, coordinates, 4584, 3684, 3819,
                                                        4044, 5, nmax);

        simdtrf::compute_hrr_geom_010y_dp_out_of_second(buffer, coordinates, 4674, 3729, 3819,
                                                        4134, 5, nmax);

        simdtrf::compute_hrr_geom_010z_dp_out_of_second(buffer, coordinates, 4764, 3774, 3819,
                                                        4224, 5, nmax);

        simdtrf::compute_hrr_geom_010x_dp_out_of_second(buffer, coordinates, 4854, 3864, 3999,
                                                        4314, 5, nmax);

        simdtrf::compute_hrr_geom_010y_dp_out_of_second(buffer, coordinates, 4944, 3909, 3999,
                                                        4404, 5, nmax);

        simdtrf::compute_hrr_geom_010z_dp_out_of_second(buffer, coordinates, 5034, 3954, 3999,
                                                        4494, 5, nmax);

        simdtrf::transform_p_inner(buffer, 5124, 4854, 6, 5, nmax);

        simdtrf::transform_d_outer(values + n * npairs, nvalues, buffer, 5124, 15, nmax);

        simdtrf::transform_p_inner(buffer, 5124, 4944, 6, 5, nmax);

        simdtrf::transform_d_outer(values + 75 * nvalues + n * npairs, nvalues, buffer, 5124, 15,
                                   nmax);

        simdtrf::transform_p_inner(buffer, 5124, 5034, 6, 5, nmax);

        simdtrf::transform_d_outer(values + 150 * nvalues + n * npairs, nvalues, buffer, 5124,
                                   15, nmax);

        simdtrf::transform_p_inner(buffer, 5124, 4584, 6, 5, nmax);

        simdtrf::transform_d_outer(values + 225 * nvalues + n * npairs, nvalues, buffer, 5124,
                                   15, nmax);

        simdtrf::transform_p_inner(buffer, 5124, 4674, 6, 5, nmax);

        simdtrf::transform_d_outer(values + 300 * nvalues + n * npairs, nvalues, buffer, 5124,
                                   15, nmax);

        simdtrf::transform_p_inner(buffer, 5124, 4764, 6, 5, nmax);

        simdtrf::transform_d_outer(values + 375 * nvalues + n * npairs, nvalues, buffer, 5124,
                                   15, nmax);
    }

    for (size_t m = 0; m < 450; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
