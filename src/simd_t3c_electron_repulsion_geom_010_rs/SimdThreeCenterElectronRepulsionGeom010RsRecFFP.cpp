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


#include "SimdThreeCenterElectronRepulsionGeom010RsRecFFP.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

#include "SimdGeometryF1.hpp"
#include "SimdGeometryG1.hpp"
#include "SimdGeometryH1.hpp"
#include "SimdGeometryI1.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdTransferDF.hpp"
#include "SimdTransferGeom010XDF.hpp"
#include "SimdTransferGeom010XDG.hpp"
#include "SimdTransferGeom010XFF.hpp"
#include "SimdTransferGeom010XPF.hpp"
#include "SimdTransferGeom010XPG.hpp"
#include "SimdTransferGeom010XPH.hpp"
#include "SimdTransferGeom010YDF.hpp"
#include "SimdTransferGeom010YDG.hpp"
#include "SimdTransferGeom010YFF.hpp"
#include "SimdTransferGeom010YPF.hpp"
#include "SimdTransferGeom010YPG.hpp"
#include "SimdTransferGeom010YPH.hpp"
#include "SimdTransferGeom010ZDF.hpp"
#include "SimdTransferGeom010ZDG.hpp"
#include "SimdTransferGeom010ZFF.hpp"
#include "SimdTransferGeom010ZPF.hpp"
#include "SimdTransferGeom010ZPG.hpp"
#include "SimdTransferGeom010ZPH.hpp"
#include "SimdTransferPF.hpp"
#include "SimdTransferPG.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_010_ffp_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_010_ffp_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 13916, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 882 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 13916, 2696, 3132, dimensions);

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

                    simdfunc::compute_t3c_erf_boys_function(buffer, coordinates, 6, 3, {1, 2, 3,
                                                            4, 5, 6, 7, 8}, ncols, fj,
                                                            i * nprim_b + j, fq, omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 15, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8}, ncols, fj, i * nprim_b + j,
                                                        fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 24, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 27, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 30, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 33, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 36, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 39, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 42, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 45, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 48, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 51, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 54, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 57, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 60, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 63, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 66, 0, 3, 7, 8,
                                                                       24, 27, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 72, 0, 3, 8, 9,
                                                                       27, 30, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 78, 0, 3, 9, 10,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 84, 0, 3, 10, 11,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 90, 0, 3, 11, 12,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 96, 0, 3, 12, 13,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 102, 0, 3, 16, 17,
                                                                       45, 48, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 108, 0, 3, 17, 18,
                                                                       48, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 114, 0, 3, 18, 19,
                                                                       51, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 120, 0, 3, 19, 20,
                                                                       54, 57, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 126, 0, 3, 20, 21,
                                                                       57, 60, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 132, 0, 3, 21, 22,
                                                                       60, 63, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 138, 0, 3, 24, 27,
                                                                       66, 72, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 148, 0, 3, 27, 30,
                                                                       72, 78, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 158, 0, 3, 30, 33,
                                                                       78, 84, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 168, 0, 3, 33, 36,
                                                                       84, 90, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 178, 0, 3, 36, 39,
                                                                       90, 96, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 188, 0, 3, 45, 48,
                                                                       102, 108, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 198, 0, 3, 48, 51,
                                                                       108, 114, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 208, 0, 3, 51, 54,
                                                                       114, 120, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 218, 0, 3, 54, 57,
                                                                       120, 126, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 228, 0, 3, 57, 60,
                                                                       126, 132, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 238, 0, 3, 66, 72,
                                                                       138, 148, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 253, 0, 3, 72, 78,
                                                                       148, 158, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 268, 0, 3, 78, 84,
                                                                       158, 168, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 283, 0, 3, 84, 90,
                                                                       168, 178, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 298, 0, 3, 102,
                                                                       108, 188, 198, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 313, 0, 3, 108,
                                                                       114, 198, 208, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 328, 0, 3, 114,
                                                                       120, 208, 218, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 343, 0, 3, 120,
                                                                       126, 218, 228, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 358, 0, 3, 138,
                                                                       148, 238, 253, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 379, 0, 3, 148,
                                                                       158, 253, 268, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 400, 0, 3, 158,
                                                                       168, 268, 283, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 421, 0, 3, 188,
                                                                       198, 298, 313, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 442, 0, 3, 198,
                                                                       208, 313, 328, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 463, 0, 3, 208,
                                                                       218, 328, 343, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 484, 0, 3, 238,
                                                                       253, 358, 379, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 512, 0, 3, 253,
                                                                       268, 379, 400, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 540, 0, 3, 298,
                                                                       313, 421, 442, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 568, 0, 3, 313,
                                                                       328, 442, 463, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 596, 0, 3, 358,
                                                                       379, 484, 512, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 632, 0, 3, 421,
                                                                       442, 540, 568, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 668, 3, 24, 66,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 686, 3, 45, 102,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 704, 3, 66, 138,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 734, 3, 102, 188,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 764, 3, 138, 238,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 809, 3, 188, 298,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 854, 3, 238, 358,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 917, 3, 298, 421,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 980, 3, 358, 484,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 1064, 3, 421, 540,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 1148, 3, 484, 596,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 1256, 3, 540, 632,
                                                                       ncols, p, q);

                    simdgeo::geom_f_x(buffer, 1364, 668, 764, 1, 3, ncols, beta);

                    simdgeo::geom_f_y(buffer, 1394, 668, 764, 1, 3, ncols, beta);

                    simdgeo::geom_f_z(buffer, 1424, 668, 764, 1, 3, ncols, beta);

                    simdgeo::geom_f_x(buffer, 1454, 686, 809, 1, 3, ncols, beta);

                    simdgeo::geom_f_y(buffer, 1484, 686, 809, 1, 3, ncols, beta);

                    simdgeo::geom_f_z(buffer, 1514, 686, 809, 1, 3, ncols, beta);

                    simdgeo::geom_g_x(buffer, 1544, 704, 854, 1, 3, ncols, beta);

                    simdgeo::geom_g_y(buffer, 1589, 704, 854, 1, 3, ncols, beta);

                    simdgeo::geom_g_z(buffer, 1634, 704, 854, 1, 3, ncols, beta);

                    simdgeo::geom_g_x(buffer, 1679, 734, 917, 1, 3, ncols, beta);

                    simdgeo::geom_g_y(buffer, 1724, 734, 917, 1, 3, ncols, beta);

                    simdgeo::geom_g_z(buffer, 1769, 734, 917, 1, 3, ncols, beta);

                    simdgeo::geom_h_x(buffer, 1814, 764, 980, 1, 3, ncols, beta);

                    simdgeo::geom_h_y(buffer, 1877, 764, 980, 1, 3, ncols, beta);

                    simdgeo::geom_h_z(buffer, 1940, 764, 980, 1, 3, ncols, beta);

                    simdgeo::geom_h_x(buffer, 2003, 809, 1064, 1, 3, ncols, beta);

                    simdgeo::geom_h_y(buffer, 2066, 809, 1064, 1, 3, ncols, beta);

                    simdgeo::geom_h_z(buffer, 2129, 809, 1064, 1, 3, ncols, beta);

                    simdgeo::geom_i_x(buffer, 2192, 854, 1148, 1, 3, ncols, beta);

                    simdgeo::geom_i_y(buffer, 2276, 854, 1148, 1, 3, ncols, beta);

                    simdgeo::geom_i_z(buffer, 2360, 854, 1148, 1, 3, ncols, beta);

                    simdgeo::geom_i_x(buffer, 2444, 917, 1256, 1, 3, ncols, beta);

                    simdgeo::geom_i_y(buffer, 2528, 917, 1256, 1, 3, ncols, beta);

                    simdgeo::geom_i_z(buffer, 2612, 917, 1256, 1, 3, ncols, beta);

                    simdfunc::contract_primitives(buffer, 2696, 1364, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2756, 1394, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2816, 1424, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2876, 704, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2936, 1454, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2996, 1484, 30, ncols);

                    simdfunc::contract_primitives(buffer, 3056, 1514, 30, ncols);

                    simdfunc::contract_primitives(buffer, 3116, 734, 30, ncols);

                    simdfunc::contract_primitives(buffer, 3176, 1544, 45, ncols);

                    simdfunc::contract_primitives(buffer, 3266, 1589, 45, ncols);

                    simdfunc::contract_primitives(buffer, 3356, 1634, 45, ncols);

                    simdfunc::contract_primitives(buffer, 3446, 764, 45, ncols);

                    simdfunc::contract_primitives(buffer, 3536, 1679, 45, ncols);

                    simdfunc::contract_primitives(buffer, 3626, 1724, 45, ncols);

                    simdfunc::contract_primitives(buffer, 3716, 1769, 45, ncols);

                    simdfunc::contract_primitives(buffer, 3806, 809, 45, ncols);

                    simdfunc::contract_primitives(buffer, 3896, 1814, 63, ncols);

                    simdfunc::contract_primitives(buffer, 4022, 1877, 63, ncols);

                    simdfunc::contract_primitives(buffer, 4148, 1940, 63, ncols);

                    simdfunc::contract_primitives(buffer, 4274, 854, 63, ncols);

                    simdfunc::contract_primitives(buffer, 4400, 2003, 63, ncols);

                    simdfunc::contract_primitives(buffer, 4526, 2066, 63, ncols);

                    simdfunc::contract_primitives(buffer, 4652, 2129, 63, ncols);

                    simdfunc::contract_primitives(buffer, 4778, 917, 63, ncols);

                    simdfunc::contract_primitives(buffer, 4904, 2192, 84, ncols);

                    simdfunc::contract_primitives(buffer, 5072, 2276, 84, ncols);

                    simdfunc::contract_primitives(buffer, 5240, 2360, 84, ncols);

                    simdfunc::contract_primitives(buffer, 5408, 2444, 84, ncols);

                    simdfunc::contract_primitives(buffer, 5576, 2528, 84, ncols);

                    simdfunc::contract_primitives(buffer, 5744, 2612, 84, ncols);
                }
            }
        }

        simdtrf::transform_p_inner(buffer, 2726, 2696, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2786, 2756, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2846, 2816, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2906, 2876, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2966, 2936, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3026, 2996, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3086, 3056, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3146, 3116, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3221, 3176, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3311, 3266, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3401, 3356, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3491, 3446, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3581, 3536, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3671, 3626, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3761, 3716, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3851, 3806, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3959, 3896, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4085, 4022, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4211, 4148, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4337, 4274, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4463, 4400, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4589, 4526, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4715, 4652, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4841, 4778, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4988, 4904, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 5156, 5072, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 5324, 5240, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 5492, 5408, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 5660, 5576, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 5828, 5744, 28, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 5912, 2726, 2906, 3221, 3, nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 6002, 2786, 2906, 3311, 3, nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 6092, 2846, 2906, 3401, 3, nmax);

        simdtrf::compute_hrr_pf(buffer, coordinates, 6182, 2906, 3491, 3, nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 6272, 2966, 3146, 3581, 3, nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 6362, 3026, 3146, 3671, 3, nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 6452, 3086, 3146, 3761, 3, nmax);

        simdtrf::compute_hrr_pf(buffer, coordinates, 6542, 3146, 3851, 3, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 6632, 3221, 3491, 3959, 3, nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 6767, 3311, 3491, 4085, 3, nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 6902, 3401, 3491, 4211, 3, nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 7037, 3491, 4337, 3, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 7172, 3581, 3851, 4463, 3, nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 7307, 3671, 3851, 4589, 3, nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 7442, 3761, 3851, 4715, 3, nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 7577, 3851, 4841, 3, nmax);

        simdtrf::compute_hrr_geom_010x_ph(buffer, coordinates, 7712, 3959, 4337, 4988, 3, nmax);

        simdtrf::compute_hrr_geom_010y_ph(buffer, coordinates, 7901, 4085, 4337, 5156, 3, nmax);

        simdtrf::compute_hrr_geom_010z_ph(buffer, coordinates, 8090, 4211, 4337, 5324, 3, nmax);

        simdtrf::compute_hrr_geom_010x_ph(buffer, coordinates, 8279, 4463, 4841, 5492, 3, nmax);

        simdtrf::compute_hrr_geom_010y_ph(buffer, coordinates, 8468, 4589, 4841, 5660, 3, nmax);

        simdtrf::compute_hrr_geom_010z_ph(buffer, coordinates, 8657, 4715, 4841, 5828, 3, nmax);

        simdtrf::compute_hrr_geom_010x_df(buffer, coordinates, 8846, 5912, 6182, 6632, 3, nmax);

        simdtrf::compute_hrr_geom_010y_df(buffer, coordinates, 9026, 6002, 6182, 6767, 3, nmax);

        simdtrf::compute_hrr_geom_010z_df(buffer, coordinates, 9206, 6092, 6182, 6902, 3, nmax);

        simdtrf::compute_hrr_df(buffer, coordinates, 9386, 6182, 7037, 3, nmax);

        simdtrf::compute_hrr_geom_010x_df(buffer, coordinates, 9566, 6272, 6542, 7172, 3, nmax);

        simdtrf::compute_hrr_geom_010y_df(buffer, coordinates, 9746, 6362, 6542, 7307, 3, nmax);

        simdtrf::compute_hrr_geom_010z_df(buffer, coordinates, 9926, 6452, 6542, 7442, 3, nmax);

        simdtrf::compute_hrr_df(buffer, coordinates, 10106, 6542, 7577, 3, nmax);

        simdtrf::compute_hrr_geom_010x_dg(buffer, coordinates, 10286, 6632, 7037, 7712, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dg(buffer, coordinates, 10556, 6767, 7037, 7901, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dg(buffer, coordinates, 10826, 6902, 7037, 8090, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_dg(buffer, coordinates, 11096, 7172, 7577, 8279, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dg(buffer, coordinates, 11366, 7307, 7577, 8468, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dg(buffer, coordinates, 11636, 7442, 7577, 8657, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_ff(buffer, coordinates, 11906, 8846, 9386, 10286, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_ff(buffer, coordinates, 12206, 9026, 9386, 10556, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_ff(buffer, coordinates, 12506, 9206, 9386, 10826, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_ff(buffer, coordinates, 12806, 9566, 10106, 11096, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_ff(buffer, coordinates, 13106, 9746, 10106, 11366, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_ff(buffer, coordinates, 13406, 9926, 10106, 11636, 3,
                                          nmax);

        simdtrf::transform_f_inner(buffer, 13706, 12806, 10, 3, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 13706, 21, nmax);

        simdtrf::transform_f_inner(buffer, 13706, 13106, 10, 3, nmax);

        simdtrf::transform_f_outer(values + 147 * nvalues + n * npairs, nvalues, buffer, 13706,
                                   21, nmax);

        simdtrf::transform_f_inner(buffer, 13706, 13406, 10, 3, nmax);

        simdtrf::transform_f_outer(values + 294 * nvalues + n * npairs, nvalues, buffer, 13706,
                                   21, nmax);

        simdtrf::transform_f_inner(buffer, 13706, 11906, 10, 3, nmax);

        simdtrf::transform_f_outer(values + 441 * nvalues + n * npairs, nvalues, buffer, 13706,
                                   21, nmax);

        simdtrf::transform_f_inner(buffer, 13706, 12206, 10, 3, nmax);

        simdtrf::transform_f_outer(values + 588 * nvalues + n * npairs, nvalues, buffer, 13706,
                                   21, nmax);

        simdtrf::transform_f_inner(buffer, 13706, 12506, 10, 3, nmax);

        simdtrf::transform_f_outer(values + 735 * nvalues + n * npairs, nvalues, buffer, 13706,
                                   21, nmax);
    }

    for (size_t m = 0; m < 882; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
