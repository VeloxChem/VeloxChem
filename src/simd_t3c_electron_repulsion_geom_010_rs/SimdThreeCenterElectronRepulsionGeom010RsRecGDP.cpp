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


#include "SimdThreeCenterElectronRepulsionGeom010RsRecGDP.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdTransferDD.hpp"
#include "SimdTransferDF.hpp"
#include "SimdTransferFD.hpp"
#include "SimdTransferGeom010XDD.hpp"
#include "SimdTransferGeom010XDF.hpp"
#include "SimdTransferGeom010XDG.hpp"
#include "SimdTransferGeom010XFD.hpp"
#include "SimdTransferGeom010XFF.hpp"
#include "SimdTransferGeom010XGD.hpp"
#include "SimdTransferGeom010XPD.hpp"
#include "SimdTransferGeom010XPF.hpp"
#include "SimdTransferGeom010XPG.hpp"
#include "SimdTransferGeom010XPH.hpp"
#include "SimdTransferGeom010YDD.hpp"
#include "SimdTransferGeom010YDF.hpp"
#include "SimdTransferGeom010YDG.hpp"
#include "SimdTransferGeom010YFD.hpp"
#include "SimdTransferGeom010YFF.hpp"
#include "SimdTransferGeom010YGD.hpp"
#include "SimdTransferGeom010YPD.hpp"
#include "SimdTransferGeom010YPF.hpp"
#include "SimdTransferGeom010YPG.hpp"
#include "SimdTransferGeom010YPH.hpp"
#include "SimdTransferGeom010ZDD.hpp"
#include "SimdTransferGeom010ZDF.hpp"
#include "SimdTransferGeom010ZDG.hpp"
#include "SimdTransferGeom010ZFD.hpp"
#include "SimdTransferGeom010ZFF.hpp"
#include "SimdTransferGeom010ZGD.hpp"
#include "SimdTransferGeom010ZPD.hpp"
#include "SimdTransferGeom010ZPF.hpp"
#include "SimdTransferGeom010ZPG.hpp"
#include "SimdTransferGeom010ZPH.hpp"
#include "SimdTransferPD.hpp"
#include "SimdTransferPF.hpp"
#include "SimdTransferPG.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_010_gdp_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_010_gdp_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 18701, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 810 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 18701, 2822, 3420, dimensions);

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

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 668, 3, 7, 24,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 677, 3, 16, 45,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 686, 3, 24, 66,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 704, 3, 45, 102,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 722, 3, 66, 138,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 752, 3, 102, 188,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 782, 3, 138, 238,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 827, 3, 188, 298,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 872, 3, 238, 358,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 935, 3, 298, 421,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 998, 3, 358, 484,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 1082, 3, 421, 540,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 1166, 3, 484, 596,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 1274, 3, 540, 632,
                                                                       ncols, p, q);

                    simdgeo::geom_d_x(buffer, 1382, 668, 722, 1, 3, ncols, beta);

                    simdgeo::geom_d_y(buffer, 1400, 668, 722, 1, 3, ncols, beta);

                    simdgeo::geom_d_z(buffer, 1418, 668, 722, 1, 3, ncols, beta);

                    simdgeo::geom_d_x(buffer, 1436, 677, 752, 1, 3, ncols, beta);

                    simdgeo::geom_d_y(buffer, 1454, 677, 752, 1, 3, ncols, beta);

                    simdgeo::geom_d_z(buffer, 1472, 677, 752, 1, 3, ncols, beta);

                    simdgeo::geom_f_x(buffer, 1490, 686, 782, 1, 3, ncols, beta);

                    simdgeo::geom_f_y(buffer, 1520, 686, 782, 1, 3, ncols, beta);

                    simdgeo::geom_f_z(buffer, 1550, 686, 782, 1, 3, ncols, beta);

                    simdgeo::geom_f_x(buffer, 1580, 704, 827, 1, 3, ncols, beta);

                    simdgeo::geom_f_y(buffer, 1610, 704, 827, 1, 3, ncols, beta);

                    simdgeo::geom_f_z(buffer, 1640, 704, 827, 1, 3, ncols, beta);

                    simdgeo::geom_g_x(buffer, 1670, 722, 872, 1, 3, ncols, beta);

                    simdgeo::geom_g_y(buffer, 1715, 722, 872, 1, 3, ncols, beta);

                    simdgeo::geom_g_z(buffer, 1760, 722, 872, 1, 3, ncols, beta);

                    simdgeo::geom_g_x(buffer, 1805, 752, 935, 1, 3, ncols, beta);

                    simdgeo::geom_g_y(buffer, 1850, 752, 935, 1, 3, ncols, beta);

                    simdgeo::geom_g_z(buffer, 1895, 752, 935, 1, 3, ncols, beta);

                    simdgeo::geom_h_x(buffer, 1940, 782, 998, 1, 3, ncols, beta);

                    simdgeo::geom_h_y(buffer, 2003, 782, 998, 1, 3, ncols, beta);

                    simdgeo::geom_h_z(buffer, 2066, 782, 998, 1, 3, ncols, beta);

                    simdgeo::geom_h_x(buffer, 2129, 827, 1082, 1, 3, ncols, beta);

                    simdgeo::geom_h_y(buffer, 2192, 827, 1082, 1, 3, ncols, beta);

                    simdgeo::geom_h_z(buffer, 2255, 827, 1082, 1, 3, ncols, beta);

                    simdgeo::geom_i_x(buffer, 2318, 872, 1166, 1, 3, ncols, beta);

                    simdgeo::geom_i_y(buffer, 2402, 872, 1166, 1, 3, ncols, beta);

                    simdgeo::geom_i_z(buffer, 2486, 872, 1166, 1, 3, ncols, beta);

                    simdgeo::geom_i_x(buffer, 2570, 935, 1274, 1, 3, ncols, beta);

                    simdgeo::geom_i_y(buffer, 2654, 935, 1274, 1, 3, ncols, beta);

                    simdgeo::geom_i_z(buffer, 2738, 935, 1274, 1, 3, ncols, beta);

                    simdfunc::contract_primitives(buffer, 2822, 1382, 18, ncols);

                    simdfunc::contract_primitives(buffer, 2858, 1400, 18, ncols);

                    simdfunc::contract_primitives(buffer, 2894, 1418, 18, ncols);

                    simdfunc::contract_primitives(buffer, 2930, 686, 18, ncols);

                    simdfunc::contract_primitives(buffer, 2966, 1436, 18, ncols);

                    simdfunc::contract_primitives(buffer, 3002, 1454, 18, ncols);

                    simdfunc::contract_primitives(buffer, 3038, 1472, 18, ncols);

                    simdfunc::contract_primitives(buffer, 3074, 704, 18, ncols);

                    simdfunc::contract_primitives(buffer, 3110, 1490, 30, ncols);

                    simdfunc::contract_primitives(buffer, 3170, 1520, 30, ncols);

                    simdfunc::contract_primitives(buffer, 3230, 1550, 30, ncols);

                    simdfunc::contract_primitives(buffer, 3290, 722, 30, ncols);

                    simdfunc::contract_primitives(buffer, 3350, 1580, 30, ncols);

                    simdfunc::contract_primitives(buffer, 3410, 1610, 30, ncols);

                    simdfunc::contract_primitives(buffer, 3470, 1640, 30, ncols);

                    simdfunc::contract_primitives(buffer, 3530, 752, 30, ncols);

                    simdfunc::contract_primitives(buffer, 3590, 1670, 45, ncols);

                    simdfunc::contract_primitives(buffer, 3680, 1715, 45, ncols);

                    simdfunc::contract_primitives(buffer, 3770, 1760, 45, ncols);

                    simdfunc::contract_primitives(buffer, 3860, 782, 45, ncols);

                    simdfunc::contract_primitives(buffer, 3950, 1805, 45, ncols);

                    simdfunc::contract_primitives(buffer, 4040, 1850, 45, ncols);

                    simdfunc::contract_primitives(buffer, 4130, 1895, 45, ncols);

                    simdfunc::contract_primitives(buffer, 4220, 827, 45, ncols);

                    simdfunc::contract_primitives(buffer, 4310, 1940, 63, ncols);

                    simdfunc::contract_primitives(buffer, 4436, 2003, 63, ncols);

                    simdfunc::contract_primitives(buffer, 4562, 2066, 63, ncols);

                    simdfunc::contract_primitives(buffer, 4688, 872, 63, ncols);

                    simdfunc::contract_primitives(buffer, 4814, 2129, 63, ncols);

                    simdfunc::contract_primitives(buffer, 4940, 2192, 63, ncols);

                    simdfunc::contract_primitives(buffer, 5066, 2255, 63, ncols);

                    simdfunc::contract_primitives(buffer, 5192, 935, 63, ncols);

                    simdfunc::contract_primitives(buffer, 5318, 2318, 84, ncols);

                    simdfunc::contract_primitives(buffer, 5486, 2402, 84, ncols);

                    simdfunc::contract_primitives(buffer, 5654, 2486, 84, ncols);

                    simdfunc::contract_primitives(buffer, 5822, 2570, 84, ncols);

                    simdfunc::contract_primitives(buffer, 5990, 2654, 84, ncols);

                    simdfunc::contract_primitives(buffer, 6158, 2738, 84, ncols);
                }
            }
        }

        simdtrf::transform_p_inner(buffer, 2840, 2822, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2876, 2858, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2912, 2894, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2948, 2930, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2984, 2966, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3020, 3002, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3056, 3038, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3092, 3074, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3140, 3110, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3200, 3170, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3260, 3230, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3320, 3290, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3380, 3350, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3440, 3410, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3500, 3470, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3560, 3530, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3635, 3590, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3725, 3680, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3815, 3770, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3905, 3860, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3995, 3950, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4085, 4040, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4175, 4130, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4265, 4220, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4373, 4310, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4499, 4436, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4625, 4562, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4751, 4688, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4877, 4814, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 5003, 4940, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 5129, 5066, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 5255, 5192, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 5402, 5318, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 5570, 5486, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 5738, 5654, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 5906, 5822, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 6074, 5990, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 6242, 6158, 28, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pd(buffer, coordinates, 6326, 2840, 2948, 3140, 3, nmax);

        simdtrf::compute_hrr_geom_010y_pd(buffer, coordinates, 6380, 2876, 2948, 3200, 3, nmax);

        simdtrf::compute_hrr_geom_010z_pd(buffer, coordinates, 6434, 2912, 2948, 3260, 3, nmax);

        simdtrf::compute_hrr_pd(buffer, coordinates, 6488, 2948, 3320, 3, nmax);

        simdtrf::compute_hrr_geom_010x_pd(buffer, coordinates, 6542, 2984, 3092, 3380, 3, nmax);

        simdtrf::compute_hrr_geom_010y_pd(buffer, coordinates, 6596, 3020, 3092, 3440, 3, nmax);

        simdtrf::compute_hrr_geom_010z_pd(buffer, coordinates, 6650, 3056, 3092, 3500, 3, nmax);

        simdtrf::compute_hrr_pd(buffer, coordinates, 6704, 3092, 3560, 3, nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 6758, 3140, 3320, 3635, 3, nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 6848, 3200, 3320, 3725, 3, nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 6938, 3260, 3320, 3815, 3, nmax);

        simdtrf::compute_hrr_pf(buffer, coordinates, 7028, 3320, 3905, 3, nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 7118, 3380, 3560, 3995, 3, nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 7208, 3440, 3560, 4085, 3, nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 7298, 3500, 3560, 4175, 3, nmax);

        simdtrf::compute_hrr_pf(buffer, coordinates, 7388, 3560, 4265, 3, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 7478, 3635, 3905, 4373, 3, nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 7613, 3725, 3905, 4499, 3, nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 7748, 3815, 3905, 4625, 3, nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 7883, 3905, 4751, 3, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 8018, 3995, 4265, 4877, 3, nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 8153, 4085, 4265, 5003, 3, nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 8288, 4175, 4265, 5129, 3, nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 8423, 4265, 5255, 3, nmax);

        simdtrf::compute_hrr_geom_010x_ph(buffer, coordinates, 8558, 4373, 4751, 5402, 3, nmax);

        simdtrf::compute_hrr_geom_010y_ph(buffer, coordinates, 8747, 4499, 4751, 5570, 3, nmax);

        simdtrf::compute_hrr_geom_010z_ph(buffer, coordinates, 8936, 4625, 4751, 5738, 3, nmax);

        simdtrf::compute_hrr_geom_010x_ph(buffer, coordinates, 9125, 4877, 5255, 5906, 3, nmax);

        simdtrf::compute_hrr_geom_010y_ph(buffer, coordinates, 9314, 5003, 5255, 6074, 3, nmax);

        simdtrf::compute_hrr_geom_010z_ph(buffer, coordinates, 9503, 5129, 5255, 6242, 3, nmax);

        simdtrf::compute_hrr_geom_010x_dd(buffer, coordinates, 9692, 6326, 6488, 6758, 3, nmax);

        simdtrf::compute_hrr_geom_010y_dd(buffer, coordinates, 9800, 6380, 6488, 6848, 3, nmax);

        simdtrf::compute_hrr_geom_010z_dd(buffer, coordinates, 9908, 6434, 6488, 6938, 3, nmax);

        simdtrf::compute_hrr_dd(buffer, coordinates, 10016, 6488, 7028, 3, nmax);

        simdtrf::compute_hrr_geom_010x_dd(buffer, coordinates, 10124, 6542, 6704, 7118, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dd(buffer, coordinates, 10232, 6596, 6704, 7208, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dd(buffer, coordinates, 10340, 6650, 6704, 7298, 3,
                                          nmax);

        simdtrf::compute_hrr_dd(buffer, coordinates, 10448, 6704, 7388, 3, nmax);

        simdtrf::compute_hrr_geom_010x_df(buffer, coordinates, 10556, 6758, 7028, 7478, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_df(buffer, coordinates, 10736, 6848, 7028, 7613, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_df(buffer, coordinates, 10916, 6938, 7028, 7748, 3,
                                          nmax);

        simdtrf::compute_hrr_df(buffer, coordinates, 11096, 7028, 7883, 3, nmax);

        simdtrf::compute_hrr_geom_010x_df(buffer, coordinates, 11276, 7118, 7388, 8018, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_df(buffer, coordinates, 11456, 7208, 7388, 8153, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_df(buffer, coordinates, 11636, 7298, 7388, 8288, 3,
                                          nmax);

        simdtrf::compute_hrr_df(buffer, coordinates, 11816, 7388, 8423, 3, nmax);

        simdtrf::compute_hrr_geom_010x_dg(buffer, coordinates, 11996, 7478, 7883, 8558, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dg(buffer, coordinates, 12266, 7613, 7883, 8747, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dg(buffer, coordinates, 12536, 7748, 7883, 8936, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_dg(buffer, coordinates, 12806, 8018, 8423, 9125, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dg(buffer, coordinates, 13076, 8153, 8423, 9314, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dg(buffer, coordinates, 13346, 8288, 8423, 9503, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_fd_out_of_second(buffer, coordinates, 13616, 9692, 10016,
                                                        10556, 3, nmax);

        simdtrf::compute_hrr_geom_010y_fd_out_of_second(buffer, coordinates, 13796, 9800, 10016,
                                                        10736, 3, nmax);

        simdtrf::compute_hrr_geom_010z_fd_out_of_second(buffer, coordinates, 13976, 9908, 10016,
                                                        10916, 3, nmax);

        simdtrf::compute_hrr_fd_out_of_second(buffer, coordinates, 14156, 10016, 11096, 3,
                                              nmax);

        simdtrf::compute_hrr_geom_010x_fd_out_of_second(buffer, coordinates, 14336, 10124, 10448,
                                                        11276, 3, nmax);

        simdtrf::compute_hrr_geom_010y_fd_out_of_second(buffer, coordinates, 14516, 10232, 10448,
                                                        11456, 3, nmax);

        simdtrf::compute_hrr_geom_010z_fd_out_of_second(buffer, coordinates, 14696, 10340, 10448,
                                                        11636, 3, nmax);

        simdtrf::compute_hrr_fd_out_of_second(buffer, coordinates, 14876, 10448, 11816, 3,
                                              nmax);

        simdtrf::compute_hrr_geom_010x_ff(buffer, coordinates, 15056, 10556, 11096, 11996, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_ff(buffer, coordinates, 15356, 10736, 11096, 12266, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_ff(buffer, coordinates, 15656, 10916, 11096, 12536, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_ff(buffer, coordinates, 15956, 11276, 11816, 12806, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_ff(buffer, coordinates, 16256, 11456, 11816, 13076, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_ff(buffer, coordinates, 16556, 11636, 11816, 13346, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_gd_out_of_second(buffer, coordinates, 16856, 13616, 14156,
                                                        15056, 3, nmax);

        simdtrf::compute_hrr_geom_010y_gd_out_of_second(buffer, coordinates, 17126, 13796, 14156,
                                                        15356, 3, nmax);

        simdtrf::compute_hrr_geom_010z_gd_out_of_second(buffer, coordinates, 17396, 13976, 14156,
                                                        15656, 3, nmax);

        simdtrf::compute_hrr_geom_010x_gd_out_of_second(buffer, coordinates, 17666, 14336, 14876,
                                                        15956, 3, nmax);

        simdtrf::compute_hrr_geom_010y_gd_out_of_second(buffer, coordinates, 17936, 14516, 14876,
                                                        16256, 3, nmax);

        simdtrf::compute_hrr_geom_010z_gd_out_of_second(buffer, coordinates, 18206, 14696, 14876,
                                                        16556, 3, nmax);

        simdtrf::transform_d_inner(buffer, 18476, 17666, 15, 3, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 18476, 15, nmax);

        simdtrf::transform_d_inner(buffer, 18476, 17936, 15, 3, nmax);

        simdtrf::transform_g_outer(values + 135 * nvalues + n * npairs, nvalues, buffer, 18476,
                                   15, nmax);

        simdtrf::transform_d_inner(buffer, 18476, 18206, 15, 3, nmax);

        simdtrf::transform_g_outer(values + 270 * nvalues + n * npairs, nvalues, buffer, 18476,
                                   15, nmax);

        simdtrf::transform_d_inner(buffer, 18476, 16856, 15, 3, nmax);

        simdtrf::transform_g_outer(values + 405 * nvalues + n * npairs, nvalues, buffer, 18476,
                                   15, nmax);

        simdtrf::transform_d_inner(buffer, 18476, 17126, 15, 3, nmax);

        simdtrf::transform_g_outer(values + 540 * nvalues + n * npairs, nvalues, buffer, 18476,
                                   15, nmax);

        simdtrf::transform_d_inner(buffer, 18476, 17396, 15, 3, nmax);

        simdtrf::transform_g_outer(values + 675 * nvalues + n * npairs, nvalues, buffer, 18476,
                                   15, nmax);
    }

    for (size_t m = 0; m < 810; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
