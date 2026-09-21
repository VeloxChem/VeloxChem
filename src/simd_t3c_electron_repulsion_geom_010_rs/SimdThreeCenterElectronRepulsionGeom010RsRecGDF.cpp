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


#include "SimdThreeCenterElectronRepulsionGeom010RsRecGDF.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSHD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSID.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
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
#include "SimdTransformF.hpp"
#include "SimdTransformG.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_010_gdf_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_010_gdf_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 58695, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1890 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 58695, 19892, 9732, dimensions);

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
                                                            4, 5, 6, 7, 8, 9, 10}, ncols, fj,
                                                            i * nprim_b + j, fq, omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 17, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9, 10}, ncols, fj,
                                                        i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 28, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 31, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 34, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 37, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 40, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 43, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 46, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 49, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 52, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 55, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 58, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 61, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 64, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 67, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 70, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 73, 0, 3, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 76, 0, 3, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 79, 0, 3, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 82, 0, 3, 7, 8,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 88, 0, 3, 8, 9,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 94, 0, 3, 9, 10,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 100, 0, 3, 10, 11,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 106, 0, 3, 11, 12,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 112, 0, 3, 12, 13,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 118, 0, 3, 13, 14,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 124, 0, 3, 14, 15,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 130, 0, 3, 18, 19,
                                                                       55, 58, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 136, 0, 3, 19, 20,
                                                                       58, 61, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 142, 0, 3, 20, 21,
                                                                       61, 64, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 148, 0, 3, 21, 22,
                                                                       64, 67, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 154, 0, 3, 22, 23,
                                                                       67, 70, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 160, 0, 3, 23, 24,
                                                                       70, 73, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 166, 0, 3, 24, 25,
                                                                       73, 76, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 172, 0, 3, 25, 26,
                                                                       76, 79, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 178, 0, 3, 28, 31,
                                                                       82, 88, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 188, 0, 3, 31, 34,
                                                                       88, 94, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 198, 0, 3, 34, 37,
                                                                       94, 100, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 208, 0, 3, 37, 40,
                                                                       100, 106, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 218, 0, 3, 40, 43,
                                                                       106, 112, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 228, 0, 3, 43, 46,
                                                                       112, 118, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 238, 0, 3, 46, 49,
                                                                       118, 124, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 248, 0, 3, 55, 58,
                                                                       130, 136, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 258, 0, 3, 58, 61,
                                                                       136, 142, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 268, 0, 3, 61, 64,
                                                                       142, 148, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 278, 0, 3, 64, 67,
                                                                       148, 154, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 288, 0, 3, 67, 70,
                                                                       154, 160, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 298, 0, 3, 70, 73,
                                                                       160, 166, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 308, 0, 3, 73, 76,
                                                                       166, 172, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 318, 0, 3, 82, 88,
                                                                       178, 188, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 333, 0, 3, 88, 94,
                                                                       188, 198, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 348, 0, 3, 94,
                                                                       100, 198, 208, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 363, 0, 3, 100,
                                                                       106, 208, 218, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 378, 0, 3, 106,
                                                                       112, 218, 228, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 393, 0, 3, 112,
                                                                       118, 228, 238, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 408, 0, 3, 130,
                                                                       136, 248, 258, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 423, 0, 3, 136,
                                                                       142, 258, 268, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 438, 0, 3, 142,
                                                                       148, 268, 278, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 453, 0, 3, 148,
                                                                       154, 278, 288, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 468, 0, 3, 154,
                                                                       160, 288, 298, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 483, 0, 3, 160,
                                                                       166, 298, 308, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 498, 0, 3, 178,
                                                                       188, 318, 333, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 519, 0, 3, 188,
                                                                       198, 333, 348, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 540, 0, 3, 198,
                                                                       208, 348, 363, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 561, 0, 3, 208,
                                                                       218, 363, 378, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 582, 0, 3, 218,
                                                                       228, 378, 393, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 603, 0, 3, 248,
                                                                       258, 408, 423, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 624, 0, 3, 258,
                                                                       268, 423, 438, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 645, 0, 3, 268,
                                                                       278, 438, 453, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 666, 0, 3, 278,
                                                                       288, 453, 468, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 687, 0, 3, 288,
                                                                       298, 468, 483, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 708, 0, 3, 318,
                                                                       333, 498, 519, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 736, 0, 3, 333,
                                                                       348, 519, 540, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 764, 0, 3, 348,
                                                                       363, 540, 561, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 792, 0, 3, 363,
                                                                       378, 561, 582, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 820, 0, 3, 408,
                                                                       423, 603, 624, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 848, 0, 3, 423,
                                                                       438, 624, 645, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 876, 0, 3, 438,
                                                                       453, 645, 666, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 904, 0, 3, 453,
                                                                       468, 666, 687, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 932, 0, 3, 498,
                                                                       519, 708, 736, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 968, 0, 3, 519,
                                                                       540, 736, 764, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1004, 0, 3, 540,
                                                                       561, 764, 792, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1040, 0, 3, 603,
                                                                       624, 820, 848, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1076, 0, 3, 624,
                                                                       645, 848, 876, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1112, 0, 3, 645,
                                                                       666, 876, 904, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1148, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1151, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1154, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1157, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1160, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1163, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1166, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1169, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1172, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1175, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1178, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1181, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1184, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1187, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1190, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1193, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1196, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1199, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1202, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1205, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1208, 3, 9, 34,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1217, 3, 10, 37,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1226, 3, 11, 40,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1235, 3, 12, 43,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1244, 3, 13, 46,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1253, 3, 14, 49,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1262, 3, 15, 52,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1271, 3, 20, 61,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1280, 3, 21, 64,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1289, 3, 22, 67,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1298, 3, 23, 70,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1307, 3, 24, 73,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1316, 3, 25, 76,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1325, 3, 26, 79,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1334, 3, 28, 82,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1352, 3, 31, 88,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1370, 3, 34, 94,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1388, 3, 37, 100,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1406, 3, 40, 106,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1424, 3, 43, 112,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1442, 3, 46, 118,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1460, 3, 49, 124,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1478, 3, 55, 130,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1496, 3, 58, 136,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1514, 3, 61, 142,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1532, 3, 64, 148,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1550, 3, 67, 154,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1568, 3, 70, 160,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1586, 3, 73, 166,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1604, 3, 76, 172,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1622, 3, 82, 178,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1652, 3, 88, 188,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1682, 3, 94, 198,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1712, 3, 100, 208,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1742, 3, 106, 218,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1772, 3, 112, 228,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1802, 3, 118, 238,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1832, 3, 130, 248,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1862, 3, 136, 258,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1892, 3, 142, 268,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1922, 3, 148, 278,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1952, 3, 154, 288,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1982, 3, 160, 298,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2012, 3, 166, 308,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2042, 3, 178, 318,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2087, 3, 188, 333,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2132, 3, 198, 348,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2177, 3, 208, 363,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2222, 3, 218, 378,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2267, 3, 228, 393,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2312, 3, 248, 408,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2357, 3, 258, 423,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2402, 3, 268, 438,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2447, 3, 278, 453,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2492, 3, 288, 468,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2537, 3, 298, 483,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2582, 3, 318, 498,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2645, 3, 333, 519,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2708, 3, 348, 540,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2771, 3, 363, 561,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2834, 3, 378, 582,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2897, 3, 408, 603,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2960, 3, 423, 624,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3023, 3, 438, 645,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3086, 3, 453, 666,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3149, 3, 468, 687,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3212, 3, 498, 708,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3296, 3, 519, 736,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3380, 3, 540, 764,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3464, 3, 561, 792,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3548, 3, 603, 820,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3632, 3, 624, 848,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3716, 3, 645, 876,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3800, 3, 666, 904,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 3884, 3, 708, 932,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 3992, 3, 736, 968,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4100, 3, 764,
                                                                       1004, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4208, 3, 820,
                                                                       1040, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4316, 3, 848,
                                                                       1076, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4424, 3, 876,
                                                                       1112, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4532, 3, 7, 8,
                                                                       1154, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4538, 3, 8, 9,
                                                                       1157, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4544, 3, 9, 10,
                                                                       1160, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4550, 3, 10, 11,
                                                                       1163, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4556, 3, 11, 12,
                                                                       1166, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4562, 3, 12, 13,
                                                                       1169, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4568, 3, 13, 14,
                                                                       1172, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4574, 3, 14, 15,
                                                                       1175, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4580, 3, 18, 19,
                                                                       1184, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4586, 3, 19, 20,
                                                                       1187, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4592, 3, 20, 21,
                                                                       1190, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4598, 3, 21, 22,
                                                                       1193, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4604, 3, 22, 23,
                                                                       1196, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4610, 3, 23, 24,
                                                                       1199, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4616, 3, 24, 25,
                                                                       1202, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4622, 3, 25, 26,
                                                                       1205, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4628, 0, 3, 4532,
                                                                       1154, 4538, 1208, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4646, 0, 3, 4538,
                                                                       1157, 4544, 1217, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4664, 0, 3, 4544,
                                                                       1160, 4550, 1226, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4682, 0, 3, 4550,
                                                                       1163, 4556, 1235, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4700, 0, 3, 4556,
                                                                       1166, 4562, 1244, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4718, 0, 3, 4562,
                                                                       1169, 4568, 1253, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4736, 0, 3, 4568,
                                                                       1172, 4574, 1262, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4754, 0, 3, 4580,
                                                                       1184, 4586, 1271, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4772, 0, 3, 4586,
                                                                       1187, 4592, 1280, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4790, 0, 3, 4592,
                                                                       1190, 4598, 1289, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4808, 0, 3, 4598,
                                                                       1193, 4604, 1298, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4826, 0, 3, 4604,
                                                                       1196, 4610, 1307, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4844, 0, 3, 4610,
                                                                       1199, 4616, 1316, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4862, 0, 3, 4616,
                                                                       1202, 4622, 1325, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4880, 0, 3, 4628,
                                                                       1208, 4646, 82, 88, 1370,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4916, 0, 3, 4646,
                                                                       1217, 4664, 88, 94, 1388,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4952, 0, 3, 4664,
                                                                       1226, 4682, 94, 100, 1406,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4988, 0, 3, 4682,
                                                                       1235, 4700, 100, 106,
                                                                       1424, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5024, 0, 3, 4700,
                                                                       1244, 4718, 106, 112,
                                                                       1442, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5060, 0, 3, 4718,
                                                                       1253, 4736, 112, 118,
                                                                       1460, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5096, 0, 3, 4754,
                                                                       1271, 4772, 130, 136,
                                                                       1514, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5132, 0, 3, 4772,
                                                                       1280, 4790, 136, 142,
                                                                       1532, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5168, 0, 3, 4790,
                                                                       1289, 4808, 142, 148,
                                                                       1550, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5204, 0, 3, 4808,
                                                                       1298, 4826, 148, 154,
                                                                       1568, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5240, 0, 3, 4826,
                                                                       1307, 4844, 154, 160,
                                                                       1586, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5276, 0, 3, 4844,
                                                                       1316, 4862, 160, 166,
                                                                       1604, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5312, 0, 3, 4880,
                                                                       1370, 4916, 178, 188,
                                                                       1682, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5372, 0, 3, 4916,
                                                                       1388, 4952, 188, 198,
                                                                       1712, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5432, 0, 3, 4952,
                                                                       1406, 4988, 198, 208,
                                                                       1742, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5492, 0, 3, 4988,
                                                                       1424, 5024, 208, 218,
                                                                       1772, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5552, 0, 3, 5024,
                                                                       1442, 5060, 218, 228,
                                                                       1802, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5612, 0, 3, 5096,
                                                                       1514, 5132, 248, 258,
                                                                       1892, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5672, 0, 3, 5132,
                                                                       1532, 5168, 258, 268,
                                                                       1922, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5732, 0, 3, 5168,
                                                                       1550, 5204, 268, 278,
                                                                       1952, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5792, 0, 3, 5204,
                                                                       1568, 5240, 278, 288,
                                                                       1982, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5852, 0, 3, 5240,
                                                                       1586, 5276, 288, 298,
                                                                       2012, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5912, 0, 3, 5312,
                                                                       1682, 5372, 318, 333,
                                                                       2132, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6002, 0, 3, 5372,
                                                                       1712, 5432, 333, 348,
                                                                       2177, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6092, 0, 3, 5432,
                                                                       1742, 5492, 348, 363,
                                                                       2222, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6182, 0, 3, 5492,
                                                                       1772, 5552, 363, 378,
                                                                       2267, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6272, 0, 3, 5612,
                                                                       1892, 5672, 408, 423,
                                                                       2402, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6362, 0, 3, 5672,
                                                                       1922, 5732, 423, 438,
                                                                       2447, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6452, 0, 3, 5732,
                                                                       1952, 5792, 438, 453,
                                                                       2492, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6542, 0, 3, 5792,
                                                                       1982, 5852, 453, 468,
                                                                       2537, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 6632, 0, 3, 5912,
                                                                       2132, 6002, 498, 519,
                                                                       2708, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 6758, 0, 3, 6002,
                                                                       2177, 6092, 519, 540,
                                                                       2771, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 6884, 0, 3, 6092,
                                                                       2222, 6182, 540, 561,
                                                                       2834, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 7010, 0, 3, 6272,
                                                                       2402, 6362, 603, 624,
                                                                       3023, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 7136, 0, 3, 6362,
                                                                       2447, 6452, 624, 645,
                                                                       3086, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 7262, 0, 3, 6452,
                                                                       2492, 6542, 645, 666,
                                                                       3149, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 7388, 0, 3, 6632,
                                                                       2708, 6758, 708, 736,
                                                                       3380, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 7556, 0, 3, 6758,
                                                                       2771, 6884, 736, 764,
                                                                       3464, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 7724, 0, 3, 7010,
                                                                       3023, 7136, 820, 848,
                                                                       3716, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 7892, 0, 3, 7136,
                                                                       3086, 7262, 848, 876,
                                                                       3800, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 8060, 0, 3, 7388,
                                                                       3380, 7556, 932, 968,
                                                                       4100, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 8276, 0, 3, 7724,
                                                                       3716, 7892, 1040, 1076,
                                                                       4424, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8492, 3, 1148,
                                                                       1151, 4532, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8502, 3, 1151,
                                                                       1154, 4538, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8512, 3, 1154,
                                                                       1157, 4544, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8522, 3, 1157,
                                                                       1160, 4550, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8532, 3, 1160,
                                                                       1163, 4556, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8542, 3, 1163,
                                                                       1166, 4562, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8552, 3, 1166,
                                                                       1169, 4568, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8562, 3, 1169,
                                                                       1172, 4574, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8572, 3, 1178,
                                                                       1181, 4580, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8582, 3, 1181,
                                                                       1184, 4586, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8592, 3, 1184,
                                                                       1187, 4592, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8602, 3, 1187,
                                                                       1190, 4598, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8612, 3, 1190,
                                                                       1193, 4604, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8622, 3, 1193,
                                                                       1196, 4610, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8632, 3, 1196,
                                                                       1199, 4616, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8642, 3, 1199,
                                                                       1202, 4622, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8652, 0, 3, 8492,
                                                                       4532, 8502, 4628, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8682, 0, 3, 8502,
                                                                       4538, 8512, 4646, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8712, 0, 3, 8512,
                                                                       4544, 8522, 4664, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8742, 0, 3, 8522,
                                                                       4550, 8532, 4682, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8772, 0, 3, 8532,
                                                                       4556, 8542, 4700, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8802, 0, 3, 8542,
                                                                       4562, 8552, 4718, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8832, 0, 3, 8552,
                                                                       4568, 8562, 4736, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8862, 0, 3, 8572,
                                                                       4580, 8582, 4754, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8892, 0, 3, 8582,
                                                                       4586, 8592, 4772, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8922, 0, 3, 8592,
                                                                       4592, 8602, 4790, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8952, 0, 3, 8602,
                                                                       4598, 8612, 4808, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8982, 0, 3, 8612,
                                                                       4604, 8622, 4826, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9012, 0, 3, 8622,
                                                                       4610, 8632, 4844, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9042, 0, 3, 8632,
                                                                       4616, 8642, 4862, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9072, 0, 3, 8652,
                                                                       4628, 8682, 1334, 1352,
                                                                       4880, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9132, 0, 3, 8682,
                                                                       4646, 8712, 1352, 1370,
                                                                       4916, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9192, 0, 3, 8712,
                                                                       4664, 8742, 1370, 1388,
                                                                       4952, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9252, 0, 3, 8742,
                                                                       4682, 8772, 1388, 1406,
                                                                       4988, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9312, 0, 3, 8772,
                                                                       4700, 8802, 1406, 1424,
                                                                       5024, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9372, 0, 3, 8802,
                                                                       4718, 8832, 1424, 1442,
                                                                       5060, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9432, 0, 3, 8862,
                                                                       4754, 8892, 1478, 1496,
                                                                       5096, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9492, 0, 3, 8892,
                                                                       4772, 8922, 1496, 1514,
                                                                       5132, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9552, 0, 3, 8922,
                                                                       4790, 8952, 1514, 1532,
                                                                       5168, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9612, 0, 3, 8952,
                                                                       4808, 8982, 1532, 1550,
                                                                       5204, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9672, 0, 3, 8982,
                                                                       4826, 9012, 1550, 1568,
                                                                       5240, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9732, 0, 3, 9012,
                                                                       4844, 9042, 1568, 1586,
                                                                       5276, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 9792, 0, 3, 9072,
                                                                       4880, 9132, 1622, 1652,
                                                                       5312, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 9892, 0, 3, 9132,
                                                                       4916, 9192, 1652, 1682,
                                                                       5372, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 9992, 0, 3, 9192,
                                                                       4952, 9252, 1682, 1712,
                                                                       5432, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 10092, 0, 3, 9252,
                                                                       4988, 9312, 1712, 1742,
                                                                       5492, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 10192, 0, 3, 9312,
                                                                       5024, 9372, 1742, 1772,
                                                                       5552, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 10292, 0, 3, 9432,
                                                                       5096, 9492, 1832, 1862,
                                                                       5612, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 10392, 0, 3, 9492,
                                                                       5132, 9552, 1862, 1892,
                                                                       5672, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 10492, 0, 3, 9552,
                                                                       5168, 9612, 1892, 1922,
                                                                       5732, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 10592, 0, 3, 9612,
                                                                       5204, 9672, 1922, 1952,
                                                                       5792, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 10692, 0, 3, 9672,
                                                                       5240, 9732, 1952, 1982,
                                                                       5852, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 10792, 0, 3, 9792,
                                                                       5312, 9892, 2042, 2087,
                                                                       5912, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 10942, 0, 3, 9892,
                                                                       5372, 9992, 2087, 2132,
                                                                       6002, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 11092, 0, 3, 9992,
                                                                       5432, 10092, 2132, 2177,
                                                                       6092, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 11242, 0, 3,
                                                                       10092, 5492, 10192, 2177,
                                                                       2222, 6182, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 11392, 0, 3,
                                                                       10292, 5612, 10392, 2312,
                                                                       2357, 6272, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 11542, 0, 3,
                                                                       10392, 5672, 10492, 2357,
                                                                       2402, 6362, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 11692, 0, 3,
                                                                       10492, 5732, 10592, 2402,
                                                                       2447, 6452, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 11842, 0, 3,
                                                                       10592, 5792, 10692, 2447,
                                                                       2492, 6542, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 11992, 0, 3,
                                                                       10792, 5912, 10942, 2582,
                                                                       2645, 6632, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 12202, 0, 3,
                                                                       10942, 6002, 11092, 2645,
                                                                       2708, 6758, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 12412, 0, 3,
                                                                       11092, 6092, 11242, 2708,
                                                                       2771, 6884, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 12622, 0, 3,
                                                                       11392, 6272, 11542, 2897,
                                                                       2960, 7010, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 12832, 0, 3,
                                                                       11542, 6362, 11692, 2960,
                                                                       3023, 7136, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 13042, 0, 3,
                                                                       11692, 6452, 11842, 3023,
                                                                       3086, 7262, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 13252, 0, 3,
                                                                       11992, 6632, 12202, 3212,
                                                                       3296, 7388, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 13532, 0, 3,
                                                                       12202, 6758, 12412, 3296,
                                                                       3380, 7556, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 13812, 0, 3,
                                                                       12622, 7010, 12832, 3548,
                                                                       3632, 7724, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 14092, 0, 3,
                                                                       12832, 7136, 13042, 3632,
                                                                       3716, 7892, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 14372, 0, 3,
                                                                       13252, 7388, 13532, 3884,
                                                                       3992, 8060, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 14732, 0, 3,
                                                                       13812, 7724, 14092, 4208,
                                                                       4316, 8276, ncols, gamma,
                                                                       p, q);

                    simdgeo::geom_d_x(buffer, 15092, 8652, 9792, 1, 10, ncols, beta);

                    simdgeo::geom_d_y(buffer, 15152, 8652, 9792, 1, 10, ncols, beta);

                    simdgeo::geom_d_z(buffer, 15212, 8652, 9792, 1, 10, ncols, beta);

                    simdgeo::geom_d_x(buffer, 15272, 8862, 10292, 1, 10, ncols, beta);

                    simdgeo::geom_d_y(buffer, 15332, 8862, 10292, 1, 10, ncols, beta);

                    simdgeo::geom_d_z(buffer, 15392, 8862, 10292, 1, 10, ncols, beta);

                    simdgeo::geom_f_x(buffer, 15452, 9072, 10792, 1, 10, ncols, beta);

                    simdgeo::geom_f_y(buffer, 15552, 9072, 10792, 1, 10, ncols, beta);

                    simdgeo::geom_f_z(buffer, 15652, 9072, 10792, 1, 10, ncols, beta);

                    simdgeo::geom_f_x(buffer, 15752, 9432, 11392, 1, 10, ncols, beta);

                    simdgeo::geom_f_y(buffer, 15852, 9432, 11392, 1, 10, ncols, beta);

                    simdgeo::geom_f_z(buffer, 15952, 9432, 11392, 1, 10, ncols, beta);

                    simdgeo::geom_g_x(buffer, 16052, 9792, 11992, 1, 10, ncols, beta);

                    simdgeo::geom_g_y(buffer, 16202, 9792, 11992, 1, 10, ncols, beta);

                    simdgeo::geom_g_z(buffer, 16352, 9792, 11992, 1, 10, ncols, beta);

                    simdgeo::geom_g_x(buffer, 16502, 10292, 12622, 1, 10, ncols, beta);

                    simdgeo::geom_g_y(buffer, 16652, 10292, 12622, 1, 10, ncols, beta);

                    simdgeo::geom_g_z(buffer, 16802, 10292, 12622, 1, 10, ncols, beta);

                    simdgeo::geom_h_x(buffer, 16952, 10792, 13252, 1, 10, ncols, beta);

                    simdgeo::geom_h_y(buffer, 17162, 10792, 13252, 1, 10, ncols, beta);

                    simdgeo::geom_h_z(buffer, 17372, 10792, 13252, 1, 10, ncols, beta);

                    simdgeo::geom_h_x(buffer, 17582, 11392, 13812, 1, 10, ncols, beta);

                    simdgeo::geom_h_y(buffer, 17792, 11392, 13812, 1, 10, ncols, beta);

                    simdgeo::geom_h_z(buffer, 18002, 11392, 13812, 1, 10, ncols, beta);

                    simdgeo::geom_i_x(buffer, 18212, 11992, 14372, 1, 10, ncols, beta);

                    simdgeo::geom_i_y(buffer, 18492, 11992, 14372, 1, 10, ncols, beta);

                    simdgeo::geom_i_z(buffer, 18772, 11992, 14372, 1, 10, ncols, beta);

                    simdgeo::geom_i_x(buffer, 19052, 12622, 14732, 1, 10, ncols, beta);

                    simdgeo::geom_i_y(buffer, 19332, 12622, 14732, 1, 10, ncols, beta);

                    simdgeo::geom_i_z(buffer, 19612, 12622, 14732, 1, 10, ncols, beta);

                    simdfunc::contract_primitives(buffer, 19892, 15092, 60, ncols);

                    simdfunc::contract_primitives(buffer, 19994, 15152, 60, ncols);

                    simdfunc::contract_primitives(buffer, 20096, 15212, 60, ncols);

                    simdfunc::contract_primitives(buffer, 20198, 9072, 60, ncols);

                    simdfunc::contract_primitives(buffer, 20300, 15272, 60, ncols);

                    simdfunc::contract_primitives(buffer, 20402, 15332, 60, ncols);

                    simdfunc::contract_primitives(buffer, 20504, 15392, 60, ncols);

                    simdfunc::contract_primitives(buffer, 20606, 9432, 60, ncols);

                    simdfunc::contract_primitives(buffer, 20708, 15452, 100, ncols);

                    simdfunc::contract_primitives(buffer, 20878, 15552, 100, ncols);

                    simdfunc::contract_primitives(buffer, 21048, 15652, 100, ncols);

                    simdfunc::contract_primitives(buffer, 21218, 9792, 100, ncols);

                    simdfunc::contract_primitives(buffer, 21388, 15752, 100, ncols);

                    simdfunc::contract_primitives(buffer, 21558, 15852, 100, ncols);

                    simdfunc::contract_primitives(buffer, 21728, 15952, 100, ncols);

                    simdfunc::contract_primitives(buffer, 21898, 10292, 100, ncols);

                    simdfunc::contract_primitives(buffer, 22068, 16052, 150, ncols);

                    simdfunc::contract_primitives(buffer, 22323, 16202, 150, ncols);

                    simdfunc::contract_primitives(buffer, 22578, 16352, 150, ncols);

                    simdfunc::contract_primitives(buffer, 22833, 10792, 150, ncols);

                    simdfunc::contract_primitives(buffer, 23088, 16502, 150, ncols);

                    simdfunc::contract_primitives(buffer, 23343, 16652, 150, ncols);

                    simdfunc::contract_primitives(buffer, 23598, 16802, 150, ncols);

                    simdfunc::contract_primitives(buffer, 23853, 11392, 150, ncols);

                    simdfunc::contract_primitives(buffer, 24108, 16952, 210, ncols);

                    simdfunc::contract_primitives(buffer, 24465, 17162, 210, ncols);

                    simdfunc::contract_primitives(buffer, 24822, 17372, 210, ncols);

                    simdfunc::contract_primitives(buffer, 25179, 11992, 210, ncols);

                    simdfunc::contract_primitives(buffer, 25536, 17582, 210, ncols);

                    simdfunc::contract_primitives(buffer, 25893, 17792, 210, ncols);

                    simdfunc::contract_primitives(buffer, 26250, 18002, 210, ncols);

                    simdfunc::contract_primitives(buffer, 26607, 12622, 210, ncols);

                    simdfunc::contract_primitives(buffer, 26964, 18212, 280, ncols);

                    simdfunc::contract_primitives(buffer, 27440, 18492, 280, ncols);

                    simdfunc::contract_primitives(buffer, 27916, 18772, 280, ncols);

                    simdfunc::contract_primitives(buffer, 28392, 19052, 280, ncols);

                    simdfunc::contract_primitives(buffer, 28868, 19332, 280, ncols);

                    simdfunc::contract_primitives(buffer, 29344, 19612, 280, ncols);
                }
            }
        }

        simdtrf::transform_f_inner(buffer, 19952, 19892, 6, 1, nmax);

        simdtrf::transform_f_inner(buffer, 20054, 19994, 6, 1, nmax);

        simdtrf::transform_f_inner(buffer, 20156, 20096, 6, 1, nmax);

        simdtrf::transform_f_inner(buffer, 20258, 20198, 6, 1, nmax);

        simdtrf::transform_f_inner(buffer, 20360, 20300, 6, 1, nmax);

        simdtrf::transform_f_inner(buffer, 20462, 20402, 6, 1, nmax);

        simdtrf::transform_f_inner(buffer, 20564, 20504, 6, 1, nmax);

        simdtrf::transform_f_inner(buffer, 20666, 20606, 6, 1, nmax);

        simdtrf::transform_f_inner(buffer, 20808, 20708, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 20978, 20878, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 21148, 21048, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 21318, 21218, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 21488, 21388, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 21658, 21558, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 21828, 21728, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 21998, 21898, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 22218, 22068, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 22473, 22323, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 22728, 22578, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 22983, 22833, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 23238, 23088, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 23493, 23343, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 23748, 23598, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 24003, 23853, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 24318, 24108, 21, 1, nmax);

        simdtrf::transform_f_inner(buffer, 24675, 24465, 21, 1, nmax);

        simdtrf::transform_f_inner(buffer, 25032, 24822, 21, 1, nmax);

        simdtrf::transform_f_inner(buffer, 25389, 25179, 21, 1, nmax);

        simdtrf::transform_f_inner(buffer, 25746, 25536, 21, 1, nmax);

        simdtrf::transform_f_inner(buffer, 26103, 25893, 21, 1, nmax);

        simdtrf::transform_f_inner(buffer, 26460, 26250, 21, 1, nmax);

        simdtrf::transform_f_inner(buffer, 26817, 26607, 21, 1, nmax);

        simdtrf::transform_f_inner(buffer, 27244, 26964, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 27720, 27440, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 28196, 27916, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 28672, 28392, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 29148, 28868, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 29624, 29344, 28, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pd(buffer, coordinates, 29820, 19952, 20258, 20808, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pd(buffer, coordinates, 29946, 20054, 20258, 20978, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pd(buffer, coordinates, 30072, 20156, 20258, 21148, 7,
                                          nmax);

        simdtrf::compute_hrr_pd(buffer, coordinates, 30198, 20258, 21318, 7, nmax);

        simdtrf::compute_hrr_geom_010x_pd(buffer, coordinates, 30324, 20360, 20666, 21488, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pd(buffer, coordinates, 30450, 20462, 20666, 21658, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pd(buffer, coordinates, 30576, 20564, 20666, 21828, 7,
                                          nmax);

        simdtrf::compute_hrr_pd(buffer, coordinates, 30702, 20666, 21998, 7, nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 30828, 20808, 21318, 22218, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 31038, 20978, 21318, 22473, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 31248, 21148, 21318, 22728, 7,
                                          nmax);

        simdtrf::compute_hrr_pf(buffer, coordinates, 31458, 21318, 22983, 7, nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 31668, 21488, 21998, 23238, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 31878, 21658, 21998, 23493, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 32088, 21828, 21998, 23748, 7,
                                          nmax);

        simdtrf::compute_hrr_pf(buffer, coordinates, 32298, 21998, 24003, 7, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 32508, 22218, 22983, 24318, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 32823, 22473, 22983, 24675, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 33138, 22728, 22983, 25032, 7,
                                          nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 33453, 22983, 25389, 7, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 33768, 23238, 24003, 25746, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 34083, 23493, 24003, 26103, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 34398, 23748, 24003, 26460, 7,
                                          nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 34713, 24003, 26817, 7, nmax);

        simdtrf::compute_hrr_geom_010x_ph(buffer, coordinates, 35028, 24318, 25389, 27244, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_ph(buffer, coordinates, 35469, 24675, 25389, 27720, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_ph(buffer, coordinates, 35910, 25032, 25389, 28196, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_ph(buffer, coordinates, 36351, 25746, 26817, 28672, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_ph(buffer, coordinates, 36792, 26103, 26817, 29148, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_ph(buffer, coordinates, 37233, 26460, 26817, 29624, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_dd(buffer, coordinates, 37674, 29820, 30198, 30828, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dd(buffer, coordinates, 37926, 29946, 30198, 31038, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dd(buffer, coordinates, 38178, 30072, 30198, 31248, 7,
                                          nmax);

        simdtrf::compute_hrr_dd(buffer, coordinates, 38430, 30198, 31458, 7, nmax);

        simdtrf::compute_hrr_geom_010x_dd(buffer, coordinates, 38682, 30324, 30702, 31668, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dd(buffer, coordinates, 38934, 30450, 30702, 31878, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dd(buffer, coordinates, 39186, 30576, 30702, 32088, 7,
                                          nmax);

        simdtrf::compute_hrr_dd(buffer, coordinates, 39438, 30702, 32298, 7, nmax);

        simdtrf::compute_hrr_geom_010x_df(buffer, coordinates, 39690, 30828, 31458, 32508, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_df(buffer, coordinates, 40110, 31038, 31458, 32823, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_df(buffer, coordinates, 40530, 31248, 31458, 33138, 7,
                                          nmax);

        simdtrf::compute_hrr_df(buffer, coordinates, 40950, 31458, 33453, 7, nmax);

        simdtrf::compute_hrr_geom_010x_df(buffer, coordinates, 41370, 31668, 32298, 33768, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_df(buffer, coordinates, 41790, 31878, 32298, 34083, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_df(buffer, coordinates, 42210, 32088, 32298, 34398, 7,
                                          nmax);

        simdtrf::compute_hrr_df(buffer, coordinates, 42630, 32298, 34713, 7, nmax);

        simdtrf::compute_hrr_geom_010x_dg(buffer, coordinates, 43050, 32508, 33453, 35028, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dg(buffer, coordinates, 43680, 32823, 33453, 35469, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dg(buffer, coordinates, 44310, 33138, 33453, 35910, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_dg(buffer, coordinates, 44940, 33768, 34713, 36351, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dg(buffer, coordinates, 45570, 34083, 34713, 36792, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dg(buffer, coordinates, 46200, 34398, 34713, 37233, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_fd_out_of_second(buffer, coordinates, 46830, 37674, 38430,
                                                        39690, 7, nmax);

        simdtrf::compute_hrr_geom_010y_fd_out_of_second(buffer, coordinates, 47250, 37926, 38430,
                                                        40110, 7, nmax);

        simdtrf::compute_hrr_geom_010z_fd_out_of_second(buffer, coordinates, 47670, 38178, 38430,
                                                        40530, 7, nmax);

        simdtrf::compute_hrr_fd_out_of_second(buffer, coordinates, 48090, 38430, 40950, 7,
                                              nmax);

        simdtrf::compute_hrr_geom_010x_fd_out_of_second(buffer, coordinates, 48510, 38682, 39438,
                                                        41370, 7, nmax);

        simdtrf::compute_hrr_geom_010y_fd_out_of_second(buffer, coordinates, 48930, 38934, 39438,
                                                        41790, 7, nmax);

        simdtrf::compute_hrr_geom_010z_fd_out_of_second(buffer, coordinates, 49350, 39186, 39438,
                                                        42210, 7, nmax);

        simdtrf::compute_hrr_fd_out_of_second(buffer, coordinates, 49770, 39438, 42630, 7,
                                              nmax);

        simdtrf::compute_hrr_geom_010x_ff(buffer, coordinates, 50190, 39690, 40950, 43050, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_ff(buffer, coordinates, 50890, 40110, 40950, 43680, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_ff(buffer, coordinates, 51590, 40530, 40950, 44310, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_ff(buffer, coordinates, 52290, 41370, 42630, 44940, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_ff(buffer, coordinates, 52990, 41790, 42630, 45570, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_ff(buffer, coordinates, 53690, 42210, 42630, 46200, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_gd_out_of_second(buffer, coordinates, 54390, 46830, 48090,
                                                        50190, 7, nmax);

        simdtrf::compute_hrr_geom_010y_gd_out_of_second(buffer, coordinates, 55020, 47250, 48090,
                                                        50890, 7, nmax);

        simdtrf::compute_hrr_geom_010z_gd_out_of_second(buffer, coordinates, 55650, 47670, 48090,
                                                        51590, 7, nmax);

        simdtrf::compute_hrr_geom_010x_gd_out_of_second(buffer, coordinates, 56280, 48510, 49770,
                                                        52290, 7, nmax);

        simdtrf::compute_hrr_geom_010y_gd_out_of_second(buffer, coordinates, 56910, 48930, 49770,
                                                        52990, 7, nmax);

        simdtrf::compute_hrr_geom_010z_gd_out_of_second(buffer, coordinates, 57540, 49350, 49770,
                                                        53690, 7, nmax);

        simdtrf::transform_d_inner(buffer, 58170, 56280, 15, 7, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 58170, 35, nmax);

        simdtrf::transform_d_inner(buffer, 58170, 56910, 15, 7, nmax);

        simdtrf::transform_g_outer(values + 315 * nvalues + n * npairs, nvalues, buffer, 58170,
                                   35, nmax);

        simdtrf::transform_d_inner(buffer, 58170, 57540, 15, 7, nmax);

        simdtrf::transform_g_outer(values + 630 * nvalues + n * npairs, nvalues, buffer, 58170,
                                   35, nmax);

        simdtrf::transform_d_inner(buffer, 58170, 54390, 15, 7, nmax);

        simdtrf::transform_g_outer(values + 945 * nvalues + n * npairs, nvalues, buffer, 58170,
                                   35, nmax);

        simdtrf::transform_d_inner(buffer, 58170, 55020, 15, 7, nmax);

        simdtrf::transform_g_outer(values + 1260 * nvalues + n * npairs, nvalues, buffer, 58170,
                                   35, nmax);

        simdtrf::transform_d_inner(buffer, 58170, 55650, 15, 7, nmax);

        simdtrf::transform_g_outer(values + 1575 * nvalues + n * npairs, nvalues, buffer, 58170,
                                   35, nmax);
    }

    for (size_t m = 0; m < 1890; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
