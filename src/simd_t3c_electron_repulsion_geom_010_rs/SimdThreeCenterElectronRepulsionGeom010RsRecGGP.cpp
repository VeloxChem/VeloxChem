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


#include "SimdThreeCenterElectronRepulsionGeom010RsRecGGP.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

#include "SimdGeometryG1.hpp"
#include "SimdGeometryH1.hpp"
#include "SimdGeometryI1.hpp"
#include "SimdGeometryK1.hpp"
#include "SimdGeometryL1.hpp"
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
#include "SimdThreeCenterElectronRepulsionVrrRecSLP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdTransferDG.hpp"
#include "SimdTransferDH.hpp"
#include "SimdTransferFG.hpp"
#include "SimdTransferGeom010XDG.hpp"
#include "SimdTransferGeom010XDH.hpp"
#include "SimdTransferGeom010XDI.hpp"
#include "SimdTransferGeom010XFG.hpp"
#include "SimdTransferGeom010XFH.hpp"
#include "SimdTransferGeom010XGG.hpp"
#include "SimdTransferGeom010XPG.hpp"
#include "SimdTransferGeom010XPH.hpp"
#include "SimdTransferGeom010XPI.hpp"
#include "SimdTransferGeom010XPK.hpp"
#include "SimdTransferGeom010YDG.hpp"
#include "SimdTransferGeom010YDH.hpp"
#include "SimdTransferGeom010YDI.hpp"
#include "SimdTransferGeom010YFG.hpp"
#include "SimdTransferGeom010YFH.hpp"
#include "SimdTransferGeom010YGG.hpp"
#include "SimdTransferGeom010YPG.hpp"
#include "SimdTransferGeom010YPH.hpp"
#include "SimdTransferGeom010YPI.hpp"
#include "SimdTransferGeom010YPK.hpp"
#include "SimdTransferGeom010ZDG.hpp"
#include "SimdTransferGeom010ZDH.hpp"
#include "SimdTransferGeom010ZDI.hpp"
#include "SimdTransferGeom010ZFG.hpp"
#include "SimdTransferGeom010ZFH.hpp"
#include "SimdTransferGeom010ZGG.hpp"
#include "SimdTransferGeom010ZPG.hpp"
#include "SimdTransferGeom010ZPH.hpp"
#include "SimdTransferGeom010ZPI.hpp"
#include "SimdTransferGeom010ZPK.hpp"
#include "SimdTransferPG.hpp"
#include "SimdTransferPH.hpp"
#include "SimdTransferPI.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_010_ggp_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_010_ggp_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 38323, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1458 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 38323, 5308, 6285, dimensions);

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

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1148, 0, 3, 708,
                                                                       736, 932, 968, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1193, 0, 3, 736,
                                                                       764, 968, 1004, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1238, 0, 3, 820,
                                                                       848, 1040, 1076, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1283, 0, 3, 848,
                                                                       876, 1076, 1112, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1328, 0, 3, 932,
                                                                       968, 1148, 1193, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1383, 0, 3, 1040,
                                                                       1076, 1238, 1283, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1438, 3, 82, 178,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1468, 3, 130, 248,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1498, 3, 178, 318,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1543, 3, 248, 408,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1588, 3, 318, 498,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1651, 3, 408, 603,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 1714, 3, 498, 708,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 1798, 3, 603, 820,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 1882, 3, 708, 932,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 1990, 3, 820,
                                                                       1040, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 2098, 3, 932,
                                                                       1148, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 2233, 3, 1040,
                                                                       1238, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 2368, 3, 1148,
                                                                       1328, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 2533, 3, 1238,
                                                                       1383, ncols, p, q);

                    simdgeo::geom_g_x(buffer, 2698, 1438, 1588, 1, 3, ncols, beta);

                    simdgeo::geom_g_y(buffer, 2743, 1438, 1588, 1, 3, ncols, beta);

                    simdgeo::geom_g_z(buffer, 2788, 1438, 1588, 1, 3, ncols, beta);

                    simdgeo::geom_g_x(buffer, 2833, 1468, 1651, 1, 3, ncols, beta);

                    simdgeo::geom_g_y(buffer, 2878, 1468, 1651, 1, 3, ncols, beta);

                    simdgeo::geom_g_z(buffer, 2923, 1468, 1651, 1, 3, ncols, beta);

                    simdgeo::geom_h_x(buffer, 2968, 1498, 1714, 1, 3, ncols, beta);

                    simdgeo::geom_h_y(buffer, 3031, 1498, 1714, 1, 3, ncols, beta);

                    simdgeo::geom_h_z(buffer, 3094, 1498, 1714, 1, 3, ncols, beta);

                    simdgeo::geom_h_x(buffer, 3157, 1543, 1798, 1, 3, ncols, beta);

                    simdgeo::geom_h_y(buffer, 3220, 1543, 1798, 1, 3, ncols, beta);

                    simdgeo::geom_h_z(buffer, 3283, 1543, 1798, 1, 3, ncols, beta);

                    simdgeo::geom_i_x(buffer, 3346, 1588, 1882, 1, 3, ncols, beta);

                    simdgeo::geom_i_y(buffer, 3430, 1588, 1882, 1, 3, ncols, beta);

                    simdgeo::geom_i_z(buffer, 3514, 1588, 1882, 1, 3, ncols, beta);

                    simdgeo::geom_i_x(buffer, 3598, 1651, 1990, 1, 3, ncols, beta);

                    simdgeo::geom_i_y(buffer, 3682, 1651, 1990, 1, 3, ncols, beta);

                    simdgeo::geom_i_z(buffer, 3766, 1651, 1990, 1, 3, ncols, beta);

                    simdgeo::geom_k_x(buffer, 3850, 1714, 2098, 1, 3, ncols, beta);

                    simdgeo::geom_k_y(buffer, 3958, 1714, 2098, 1, 3, ncols, beta);

                    simdgeo::geom_k_z(buffer, 4066, 1714, 2098, 1, 3, ncols, beta);

                    simdgeo::geom_k_x(buffer, 4174, 1798, 2233, 1, 3, ncols, beta);

                    simdgeo::geom_k_y(buffer, 4282, 1798, 2233, 1, 3, ncols, beta);

                    simdgeo::geom_k_z(buffer, 4390, 1798, 2233, 1, 3, ncols, beta);

                    simdgeo::geom_l_x(buffer, 4498, 1882, 2368, 1, 3, ncols, beta);

                    simdgeo::geom_l_y(buffer, 4633, 1882, 2368, 1, 3, ncols, beta);

                    simdgeo::geom_l_z(buffer, 4768, 1882, 2368, 1, 3, ncols, beta);

                    simdgeo::geom_l_x(buffer, 4903, 1990, 2533, 1, 3, ncols, beta);

                    simdgeo::geom_l_y(buffer, 5038, 1990, 2533, 1, 3, ncols, beta);

                    simdgeo::geom_l_z(buffer, 5173, 1990, 2533, 1, 3, ncols, beta);

                    simdfunc::contract_primitives(buffer, 5308, 2698, 45, ncols);

                    simdfunc::contract_primitives(buffer, 5398, 2743, 45, ncols);

                    simdfunc::contract_primitives(buffer, 5488, 2788, 45, ncols);

                    simdfunc::contract_primitives(buffer, 5578, 1498, 45, ncols);

                    simdfunc::contract_primitives(buffer, 5668, 2833, 45, ncols);

                    simdfunc::contract_primitives(buffer, 5758, 2878, 45, ncols);

                    simdfunc::contract_primitives(buffer, 5848, 2923, 45, ncols);

                    simdfunc::contract_primitives(buffer, 5938, 1543, 45, ncols);

                    simdfunc::contract_primitives(buffer, 6028, 2968, 63, ncols);

                    simdfunc::contract_primitives(buffer, 6154, 3031, 63, ncols);

                    simdfunc::contract_primitives(buffer, 6280, 3094, 63, ncols);

                    simdfunc::contract_primitives(buffer, 6406, 1588, 63, ncols);

                    simdfunc::contract_primitives(buffer, 6532, 3157, 63, ncols);

                    simdfunc::contract_primitives(buffer, 6658, 3220, 63, ncols);

                    simdfunc::contract_primitives(buffer, 6784, 3283, 63, ncols);

                    simdfunc::contract_primitives(buffer, 6910, 1651, 63, ncols);

                    simdfunc::contract_primitives(buffer, 7036, 3346, 84, ncols);

                    simdfunc::contract_primitives(buffer, 7204, 3430, 84, ncols);

                    simdfunc::contract_primitives(buffer, 7372, 3514, 84, ncols);

                    simdfunc::contract_primitives(buffer, 7540, 1714, 84, ncols);

                    simdfunc::contract_primitives(buffer, 7708, 3598, 84, ncols);

                    simdfunc::contract_primitives(buffer, 7876, 3682, 84, ncols);

                    simdfunc::contract_primitives(buffer, 8044, 3766, 84, ncols);

                    simdfunc::contract_primitives(buffer, 8212, 1798, 84, ncols);

                    simdfunc::contract_primitives(buffer, 8380, 3850, 108, ncols);

                    simdfunc::contract_primitives(buffer, 8596, 3958, 108, ncols);

                    simdfunc::contract_primitives(buffer, 8812, 4066, 108, ncols);

                    simdfunc::contract_primitives(buffer, 9028, 1882, 108, ncols);

                    simdfunc::contract_primitives(buffer, 9244, 4174, 108, ncols);

                    simdfunc::contract_primitives(buffer, 9460, 4282, 108, ncols);

                    simdfunc::contract_primitives(buffer, 9676, 4390, 108, ncols);

                    simdfunc::contract_primitives(buffer, 9892, 1990, 108, ncols);

                    simdfunc::contract_primitives(buffer, 10108, 4498, 135, ncols);

                    simdfunc::contract_primitives(buffer, 10378, 4633, 135, ncols);

                    simdfunc::contract_primitives(buffer, 10648, 4768, 135, ncols);

                    simdfunc::contract_primitives(buffer, 10918, 4903, 135, ncols);

                    simdfunc::contract_primitives(buffer, 11188, 5038, 135, ncols);

                    simdfunc::contract_primitives(buffer, 11458, 5173, 135, ncols);
                }
            }
        }

        simdtrf::transform_p_inner(buffer, 5353, 5308, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 5443, 5398, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 5533, 5488, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 5623, 5578, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 5713, 5668, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 5803, 5758, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 5893, 5848, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 5983, 5938, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 6091, 6028, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 6217, 6154, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 6343, 6280, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 6469, 6406, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 6595, 6532, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 6721, 6658, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 6847, 6784, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 6973, 6910, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 7120, 7036, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 7288, 7204, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 7456, 7372, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 7624, 7540, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 7792, 7708, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 7960, 7876, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 8128, 8044, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 8296, 8212, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 8488, 8380, 36, 1, nmax);

        simdtrf::transform_p_inner(buffer, 8704, 8596, 36, 1, nmax);

        simdtrf::transform_p_inner(buffer, 8920, 8812, 36, 1, nmax);

        simdtrf::transform_p_inner(buffer, 9136, 9028, 36, 1, nmax);

        simdtrf::transform_p_inner(buffer, 9352, 9244, 36, 1, nmax);

        simdtrf::transform_p_inner(buffer, 9568, 9460, 36, 1, nmax);

        simdtrf::transform_p_inner(buffer, 9784, 9676, 36, 1, nmax);

        simdtrf::transform_p_inner(buffer, 10000, 9892, 36, 1, nmax);

        simdtrf::transform_p_inner(buffer, 10243, 10108, 45, 1, nmax);

        simdtrf::transform_p_inner(buffer, 10513, 10378, 45, 1, nmax);

        simdtrf::transform_p_inner(buffer, 10783, 10648, 45, 1, nmax);

        simdtrf::transform_p_inner(buffer, 11053, 10918, 45, 1, nmax);

        simdtrf::transform_p_inner(buffer, 11323, 11188, 45, 1, nmax);

        simdtrf::transform_p_inner(buffer, 11593, 11458, 45, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 11728, 5353, 5623, 6091, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 11863, 5443, 5623, 6217, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 11998, 5533, 5623, 6343, 3,
                                          nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 12133, 5623, 6469, 3, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 12268, 5713, 5983, 6595, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 12403, 5803, 5983, 6721, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 12538, 5893, 5983, 6847, 3,
                                          nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 12673, 5983, 6973, 3, nmax);

        simdtrf::compute_hrr_geom_010x_ph(buffer, coordinates, 12808, 6091, 6469, 7120, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_ph(buffer, coordinates, 12997, 6217, 6469, 7288, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_ph(buffer, coordinates, 13186, 6343, 6469, 7456, 3,
                                          nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 13375, 6469, 7624, 3, nmax);

        simdtrf::compute_hrr_geom_010x_ph(buffer, coordinates, 13564, 6595, 6973, 7792, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_ph(buffer, coordinates, 13753, 6721, 6973, 7960, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_ph(buffer, coordinates, 13942, 6847, 6973, 8128, 3,
                                          nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 14131, 6973, 8296, 3, nmax);

        simdtrf::compute_hrr_geom_010x_pi(buffer, coordinates, 14320, 7120, 7624, 8488, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pi(buffer, coordinates, 14572, 7288, 7624, 8704, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pi(buffer, coordinates, 14824, 7456, 7624, 8920, 3,
                                          nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 15076, 7624, 9136, 3, nmax);

        simdtrf::compute_hrr_geom_010x_pi(buffer, coordinates, 15328, 7792, 8296, 9352, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pi(buffer, coordinates, 15580, 7960, 8296, 9568, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pi(buffer, coordinates, 15832, 8128, 8296, 9784, 3,
                                          nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 16084, 8296, 10000, 3, nmax);

        simdtrf::compute_hrr_geom_010x_pk(buffer, coordinates, 16336, 8488, 9136, 10243, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pk(buffer, coordinates, 16660, 8704, 9136, 10513, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pk(buffer, coordinates, 16984, 8920, 9136, 10783, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_pk(buffer, coordinates, 17308, 9352, 10000, 11053, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pk(buffer, coordinates, 17632, 9568, 10000, 11323, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pk(buffer, coordinates, 17956, 9784, 10000, 11593, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_dg(buffer, coordinates, 18280, 11728, 12133, 12808, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dg(buffer, coordinates, 18550, 11863, 12133, 12997, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dg(buffer, coordinates, 18820, 11998, 12133, 13186, 3,
                                          nmax);

        simdtrf::compute_hrr_dg(buffer, coordinates, 19090, 12133, 13375, 3, nmax);

        simdtrf::compute_hrr_geom_010x_dg(buffer, coordinates, 19360, 12268, 12673, 13564, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dg(buffer, coordinates, 19630, 12403, 12673, 13753, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dg(buffer, coordinates, 19900, 12538, 12673, 13942, 3,
                                          nmax);

        simdtrf::compute_hrr_dg(buffer, coordinates, 20170, 12673, 14131, 3, nmax);

        simdtrf::compute_hrr_geom_010x_dh(buffer, coordinates, 20440, 12808, 13375, 14320, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dh(buffer, coordinates, 20818, 12997, 13375, 14572, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dh(buffer, coordinates, 21196, 13186, 13375, 14824, 3,
                                          nmax);

        simdtrf::compute_hrr_dh(buffer, coordinates, 21574, 13375, 15076, 3, nmax);

        simdtrf::compute_hrr_geom_010x_dh(buffer, coordinates, 21952, 13564, 14131, 15328, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dh(buffer, coordinates, 22330, 13753, 14131, 15580, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dh(buffer, coordinates, 22708, 13942, 14131, 15832, 3,
                                          nmax);

        simdtrf::compute_hrr_dh(buffer, coordinates, 23086, 14131, 16084, 3, nmax);

        simdtrf::compute_hrr_geom_010x_di(buffer, coordinates, 23464, 14320, 15076, 16336, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_di(buffer, coordinates, 23968, 14572, 15076, 16660, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_di(buffer, coordinates, 24472, 14824, 15076, 16984, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_di(buffer, coordinates, 24976, 15328, 16084, 17308, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_di(buffer, coordinates, 25480, 15580, 16084, 17632, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_di(buffer, coordinates, 25984, 15832, 16084, 17956, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_fg(buffer, coordinates, 26488, 18280, 19090, 20440, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_fg(buffer, coordinates, 26938, 18550, 19090, 20818, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_fg(buffer, coordinates, 27388, 18820, 19090, 21196, 3,
                                          nmax);

        simdtrf::compute_hrr_fg(buffer, coordinates, 27838, 19090, 21574, 3, nmax);

        simdtrf::compute_hrr_geom_010x_fg(buffer, coordinates, 28288, 19360, 20170, 21952, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_fg(buffer, coordinates, 28738, 19630, 20170, 22330, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_fg(buffer, coordinates, 29188, 19900, 20170, 22708, 3,
                                          nmax);

        simdtrf::compute_hrr_fg(buffer, coordinates, 29638, 20170, 23086, 3, nmax);

        simdtrf::compute_hrr_geom_010x_fh(buffer, coordinates, 30088, 20440, 21574, 23464, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_fh(buffer, coordinates, 30718, 20818, 21574, 23968, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_fh(buffer, coordinates, 31348, 21196, 21574, 24472, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_fh(buffer, coordinates, 31978, 21952, 23086, 24976, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_fh(buffer, coordinates, 32608, 22330, 23086, 25480, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_fh(buffer, coordinates, 33238, 22708, 23086, 25984, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_gg(buffer, coordinates, 33868, 26488, 27838, 30088, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_gg(buffer, coordinates, 34543, 26938, 27838, 30718, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_gg(buffer, coordinates, 35218, 27388, 27838, 31348, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_gg(buffer, coordinates, 35893, 28288, 29638, 31978, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_gg(buffer, coordinates, 36568, 28738, 29638, 32608, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_gg(buffer, coordinates, 37243, 29188, 29638, 33238, 3,
                                          nmax);

        simdtrf::transform_g_inner(buffer, 37918, 35893, 15, 3, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 37918, 27, nmax);

        simdtrf::transform_g_inner(buffer, 37918, 36568, 15, 3, nmax);

        simdtrf::transform_g_outer(values + 243 * nvalues + n * npairs, nvalues, buffer, 37918,
                                   27, nmax);

        simdtrf::transform_g_inner(buffer, 37918, 37243, 15, 3, nmax);

        simdtrf::transform_g_outer(values + 486 * nvalues + n * npairs, nvalues, buffer, 37918,
                                   27, nmax);

        simdtrf::transform_g_inner(buffer, 37918, 33868, 15, 3, nmax);

        simdtrf::transform_g_outer(values + 729 * nvalues + n * npairs, nvalues, buffer, 37918,
                                   27, nmax);

        simdtrf::transform_g_inner(buffer, 37918, 34543, 15, 3, nmax);

        simdtrf::transform_g_outer(values + 972 * nvalues + n * npairs, nvalues, buffer, 37918,
                                   27, nmax);

        simdtrf::transform_g_inner(buffer, 37918, 35218, 15, 3, nmax);

        simdtrf::transform_g_outer(values + 1215 * nvalues + n * npairs, nvalues, buffer, 37918,
                                   27, nmax);
    }

    for (size_t m = 0; m < 1458; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
