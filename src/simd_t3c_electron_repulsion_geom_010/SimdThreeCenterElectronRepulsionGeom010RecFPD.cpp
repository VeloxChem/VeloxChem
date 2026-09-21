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


#include "SimdThreeCenterElectronRepulsionGeom010RecFPD.hpp"

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
#include "SimdGeometryG1.hpp"
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
#include "SimdThreeCenterElectronRepulsionVrrRecSHD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferDP.hpp"
#include "SimdTransferGeom010XDD.hpp"
#include "SimdTransferGeom010XDP.hpp"
#include "SimdTransferGeom010XFP.hpp"
#include "SimdTransferGeom010XPD.hpp"
#include "SimdTransferGeom010XPF.hpp"
#include "SimdTransferGeom010XPP.hpp"
#include "SimdTransferGeom010YDD.hpp"
#include "SimdTransferGeom010YDP.hpp"
#include "SimdTransferGeom010YFP.hpp"
#include "SimdTransferGeom010YPD.hpp"
#include "SimdTransferGeom010YPF.hpp"
#include "SimdTransferGeom010YPP.hpp"
#include "SimdTransferGeom010ZDD.hpp"
#include "SimdTransferGeom010ZDP.hpp"
#include "SimdTransferGeom010ZFP.hpp"
#include "SimdTransferGeom010ZPD.hpp"
#include "SimdTransferGeom010ZPF.hpp"
#include "SimdTransferGeom010ZPP.hpp"
#include "SimdTransferPD.hpp"
#include "SimdTransferPP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_geom_010_fpd_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_geom_010_fpd_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 5812, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 315 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 5812, 1991, 1256, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 6, 3, 7, ncols,
                                                             fj, i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 15, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 18, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 21, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 24, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 27, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 30, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 33, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 36, 0, 3, 7, 8,
                                                                       15, 18, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 42, 0, 3, 8, 9,
                                                                       18, 21, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 48, 0, 3, 9, 10,
                                                                       21, 24, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 54, 0, 3, 10, 11,
                                                                       24, 27, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 60, 0, 3, 11, 12,
                                                                       27, 30, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 66, 0, 3, 12, 13,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 72, 0, 3, 15, 18,
                                                                       36, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 82, 0, 3, 18, 21,
                                                                       42, 48, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 92, 0, 3, 21, 24,
                                                                       48, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 102, 0, 3, 24, 27,
                                                                       54, 60, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 112, 0, 3, 27, 30,
                                                                       60, 66, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 122, 0, 3, 36, 42,
                                                                       72, 82, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 137, 0, 3, 42, 48,
                                                                       82, 92, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 152, 0, 3, 48, 54,
                                                                       92, 102, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 167, 0, 3, 54, 60,
                                                                       102, 112, ncols, gamma, p,
                                                                       q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 182, 0, 3, 72, 82,
                                                                       122, 137, ncols, gamma, p,
                                                                       q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 203, 0, 3, 82, 92,
                                                                       137, 152, ncols, gamma, p,
                                                                       q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 224, 0, 3, 92,
                                                                       102, 152, 167, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 245, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 248, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 251, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 254, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 257, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 260, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 263, 3, 9, 21,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 272, 3, 10, 24,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 281, 3, 11, 27,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 290, 3, 12, 30,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 299, 3, 13, 33,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 308, 3, 21, 48,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 326, 3, 24, 54,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 344, 3, 27, 60,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 362, 3, 30, 66,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 380, 3, 48, 92,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 410, 3, 54, 102,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 440, 3, 60, 112,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 470, 3, 92, 152,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 515, 3, 102, 167,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 560, 3, 152, 224,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 623, 3, 7, 8, 245,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 629, 3, 8, 9, 248,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 635, 3, 9, 10,
                                                                       251, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 641, 3, 10, 11,
                                                                       254, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 647, 3, 11, 12,
                                                                       257, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 653, 3, 12, 13,
                                                                       260, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 659, 0, 3, 623,
                                                                       245, 629, 263, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 677, 0, 3, 629,
                                                                       248, 635, 272, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 695, 0, 3, 635,
                                                                       251, 641, 281, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 713, 0, 3, 641,
                                                                       254, 647, 290, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 731, 0, 3, 647,
                                                                       257, 653, 299, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 749, 0, 3, 659,
                                                                       263, 677, 36, 42, 308,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 785, 0, 3, 677,
                                                                       272, 695, 42, 48, 326,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 821, 0, 3, 695,
                                                                       281, 713, 48, 54, 344,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 857, 0, 3, 713,
                                                                       290, 731, 54, 60, 362,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 893, 0, 3, 749,
                                                                       308, 785, 72, 82, 380,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 953, 0, 3, 785,
                                                                       326, 821, 82, 92, 410,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 1013, 0, 3, 821,
                                                                       344, 857, 92, 102, 440,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 1073, 0, 3, 893,
                                                                       380, 953, 122, 137, 470,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 1163, 0, 3, 953,
                                                                       410, 1013, 137, 152, 515,
                                                                       ncols, gamma, p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 1253, 0, 3, 1073,
                                                                       470, 1163, 182, 203, 560,
                                                                       ncols, gamma, p, q);

                    simdgeo::geom_p_x(buffer, 1379, 623, 749, 1, 6, ncols, beta);

                    simdgeo::geom_p_y(buffer, 1397, 623, 749, 1, 6, ncols, beta);

                    simdgeo::geom_p_z(buffer, 1415, 623, 749, 1, 6, ncols, beta);

                    simdgeo::geom_d_x(buffer, 1433, 659, 893, 1, 6, ncols, beta);

                    simdgeo::geom_d_y(buffer, 1469, 659, 893, 1, 6, ncols, beta);

                    simdgeo::geom_d_z(buffer, 1505, 659, 893, 1, 6, ncols, beta);

                    simdgeo::geom_f_x(buffer, 1541, 749, 1073, 1, 6, ncols, beta);

                    simdgeo::geom_f_y(buffer, 1601, 749, 1073, 1, 6, ncols, beta);

                    simdgeo::geom_f_z(buffer, 1661, 749, 1073, 1, 6, ncols, beta);

                    simdgeo::geom_g_x(buffer, 1721, 893, 1253, 1, 6, ncols, beta);

                    simdgeo::geom_g_y(buffer, 1811, 893, 1253, 1, 6, ncols, beta);

                    simdgeo::geom_g_z(buffer, 1901, 893, 1253, 1, 6, ncols, beta);

                    simdfunc::contract_primitives(buffer, 1991, 1379, 18, ncols);

                    simdfunc::contract_primitives(buffer, 2024, 1397, 18, ncols);

                    simdfunc::contract_primitives(buffer, 2057, 1415, 18, ncols);

                    simdfunc::contract_primitives(buffer, 2090, 659, 18, ncols);

                    simdfunc::contract_primitives(buffer, 2123, 1433, 36, ncols);

                    simdfunc::contract_primitives(buffer, 2189, 1469, 36, ncols);

                    simdfunc::contract_primitives(buffer, 2255, 1505, 36, ncols);

                    simdfunc::contract_primitives(buffer, 2321, 749, 36, ncols);

                    simdfunc::contract_primitives(buffer, 2387, 1541, 60, ncols);

                    simdfunc::contract_primitives(buffer, 2497, 1601, 60, ncols);

                    simdfunc::contract_primitives(buffer, 2607, 1661, 60, ncols);

                    simdfunc::contract_primitives(buffer, 2717, 893, 60, ncols);

                    simdfunc::contract_primitives(buffer, 2827, 1721, 90, ncols);

                    simdfunc::contract_primitives(buffer, 2992, 1811, 90, ncols);

                    simdfunc::contract_primitives(buffer, 3157, 1901, 90, ncols);
                }
            }
        }

        simdtrf::transform_d_inner(buffer, 2009, 1991, 3, 1, nmax);

        simdtrf::transform_d_inner(buffer, 2042, 2024, 3, 1, nmax);

        simdtrf::transform_d_inner(buffer, 2075, 2057, 3, 1, nmax);

        simdtrf::transform_d_inner(buffer, 2108, 2090, 3, 1, nmax);

        simdtrf::transform_d_inner(buffer, 2159, 2123, 6, 1, nmax);

        simdtrf::transform_d_inner(buffer, 2225, 2189, 6, 1, nmax);

        simdtrf::transform_d_inner(buffer, 2291, 2255, 6, 1, nmax);

        simdtrf::transform_d_inner(buffer, 2357, 2321, 6, 1, nmax);

        simdtrf::transform_d_inner(buffer, 2447, 2387, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 2557, 2497, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 2667, 2607, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 2777, 2717, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 2917, 2827, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 3082, 2992, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 3247, 3157, 15, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pp(buffer, coordinates, 3322, 2009, 2108, 2159, 5, nmax);

        simdtrf::compute_hrr_geom_010y_pp(buffer, coordinates, 3367, 2042, 2108, 2225, 5, nmax);

        simdtrf::compute_hrr_geom_010z_pp(buffer, coordinates, 3412, 2075, 2108, 2291, 5, nmax);

        simdtrf::compute_hrr_pp(buffer, coordinates, 3457, 2108, 2357, 5, nmax);

        simdtrf::compute_hrr_geom_010x_pd(buffer, coordinates, 3502, 2159, 2357, 2447, 5, nmax);

        simdtrf::compute_hrr_geom_010y_pd(buffer, coordinates, 3592, 2225, 2357, 2557, 5, nmax);

        simdtrf::compute_hrr_geom_010z_pd(buffer, coordinates, 3682, 2291, 2357, 2667, 5, nmax);

        simdtrf::compute_hrr_pd(buffer, coordinates, 3772, 2357, 2777, 5, nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 3862, 2447, 2777, 2917, 5, nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 4012, 2557, 2777, 3082, 5, nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 4162, 2667, 2777, 3247, 5, nmax);

        simdtrf::compute_hrr_geom_010x_dp_out_of_second(buffer, coordinates, 4312, 3322, 3457,
                                                        3502, 5, nmax);

        simdtrf::compute_hrr_geom_010y_dp_out_of_second(buffer, coordinates, 4402, 3367, 3457,
                                                        3592, 5, nmax);

        simdtrf::compute_hrr_geom_010z_dp_out_of_second(buffer, coordinates, 4492, 3412, 3457,
                                                        3682, 5, nmax);

        simdtrf::compute_hrr_dp_out_of_second(buffer, coordinates, 4582, 3457, 3772, 5, nmax);

        simdtrf::compute_hrr_geom_010x_dd(buffer, coordinates, 4672, 3502, 3772, 3862, 5, nmax);

        simdtrf::compute_hrr_geom_010y_dd(buffer, coordinates, 4852, 3592, 3772, 4012, 5, nmax);

        simdtrf::compute_hrr_geom_010z_dd(buffer, coordinates, 5032, 3682, 3772, 4162, 5, nmax);

        simdtrf::compute_hrr_geom_010x_fp_out_of_second(buffer, coordinates, 5212, 4312, 4582,
                                                        4672, 5, nmax);

        simdtrf::compute_hrr_geom_010y_fp_out_of_second(buffer, coordinates, 5362, 4402, 4582,
                                                        4852, 5, nmax);

        simdtrf::compute_hrr_geom_010z_fp_out_of_second(buffer, coordinates, 5512, 4492, 4582,
                                                        5032, 5, nmax);

        simdtrf::transform_p_inner(buffer, 5662, 5212, 10, 5, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 5662, 15, nmax);

        simdtrf::transform_p_inner(buffer, 5662, 5362, 10, 5, nmax);

        simdtrf::transform_f_outer(values + 105 * nvalues + n * npairs, nvalues, buffer, 5662,
                                   15, nmax);

        simdtrf::transform_p_inner(buffer, 5662, 5512, 10, 5, nmax);

        simdtrf::transform_f_outer(values + 210 * nvalues + n * npairs, nvalues, buffer, 5662,
                                   15, nmax);
    }

    for (size_t m = 0; m < 315; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
