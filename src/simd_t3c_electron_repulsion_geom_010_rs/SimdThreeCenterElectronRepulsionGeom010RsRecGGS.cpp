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


#include "SimdThreeCenterElectronRepulsionGeom010RsRecGGS.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLS.hpp"
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
#include "SimdTransformS.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_010_ggs_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_010_ggs_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 13313, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 486 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 13313, 2308, 2095, dimensions);

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

                    simdfunc::compute_full_t3c_erf_boys_function(buffer, coordinates, 6, 3, 9,
                                                                 ncols, fj, i * nprim_b + j, fq,
                                                                 omega);

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 17, 3, 9,
                                                             ncols, fj, i * nprim_b + j, fq);

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

                    simdgeo::geom_g_x(buffer, 1438, 178, 498, 1, 1, ncols, beta);

                    simdgeo::geom_g_y(buffer, 1453, 178, 498, 1, 1, ncols, beta);

                    simdgeo::geom_g_z(buffer, 1468, 178, 498, 1, 1, ncols, beta);

                    simdgeo::geom_g_x(buffer, 1483, 248, 603, 1, 1, ncols, beta);

                    simdgeo::geom_g_y(buffer, 1498, 248, 603, 1, 1, ncols, beta);

                    simdgeo::geom_g_z(buffer, 1513, 248, 603, 1, 1, ncols, beta);

                    simdgeo::geom_h_x(buffer, 1528, 318, 708, 1, 1, ncols, beta);

                    simdgeo::geom_h_y(buffer, 1549, 318, 708, 1, 1, ncols, beta);

                    simdgeo::geom_h_z(buffer, 1570, 318, 708, 1, 1, ncols, beta);

                    simdgeo::geom_h_x(buffer, 1591, 408, 820, 1, 1, ncols, beta);

                    simdgeo::geom_h_y(buffer, 1612, 408, 820, 1, 1, ncols, beta);

                    simdgeo::geom_h_z(buffer, 1633, 408, 820, 1, 1, ncols, beta);

                    simdgeo::geom_i_x(buffer, 1654, 498, 932, 1, 1, ncols, beta);

                    simdgeo::geom_i_y(buffer, 1682, 498, 932, 1, 1, ncols, beta);

                    simdgeo::geom_i_z(buffer, 1710, 498, 932, 1, 1, ncols, beta);

                    simdgeo::geom_i_x(buffer, 1738, 603, 1040, 1, 1, ncols, beta);

                    simdgeo::geom_i_y(buffer, 1766, 603, 1040, 1, 1, ncols, beta);

                    simdgeo::geom_i_z(buffer, 1794, 603, 1040, 1, 1, ncols, beta);

                    simdgeo::geom_k_x(buffer, 1822, 708, 1148, 1, 1, ncols, beta);

                    simdgeo::geom_k_y(buffer, 1858, 708, 1148, 1, 1, ncols, beta);

                    simdgeo::geom_k_z(buffer, 1894, 708, 1148, 1, 1, ncols, beta);

                    simdgeo::geom_k_x(buffer, 1930, 820, 1238, 1, 1, ncols, beta);

                    simdgeo::geom_k_y(buffer, 1966, 820, 1238, 1, 1, ncols, beta);

                    simdgeo::geom_k_z(buffer, 2002, 820, 1238, 1, 1, ncols, beta);

                    simdgeo::geom_l_x(buffer, 2038, 932, 1328, 1, 1, ncols, beta);

                    simdgeo::geom_l_y(buffer, 2083, 932, 1328, 1, 1, ncols, beta);

                    simdgeo::geom_l_z(buffer, 2128, 932, 1328, 1, 1, ncols, beta);

                    simdgeo::geom_l_x(buffer, 2173, 1040, 1383, 1, 1, ncols, beta);

                    simdgeo::geom_l_y(buffer, 2218, 1040, 1383, 1, 1, ncols, beta);

                    simdgeo::geom_l_z(buffer, 2263, 1040, 1383, 1, 1, ncols, beta);

                    simdfunc::contract_primitives(buffer, 2308, 1438, 15, ncols);

                    simdfunc::contract_primitives(buffer, 2338, 1453, 15, ncols);

                    simdfunc::contract_primitives(buffer, 2368, 1468, 15, ncols);

                    simdfunc::contract_primitives(buffer, 2398, 318, 15, ncols);

                    simdfunc::contract_primitives(buffer, 2428, 1483, 15, ncols);

                    simdfunc::contract_primitives(buffer, 2458, 1498, 15, ncols);

                    simdfunc::contract_primitives(buffer, 2488, 1513, 15, ncols);

                    simdfunc::contract_primitives(buffer, 2518, 408, 15, ncols);

                    simdfunc::contract_primitives(buffer, 2548, 1528, 21, ncols);

                    simdfunc::contract_primitives(buffer, 2590, 1549, 21, ncols);

                    simdfunc::contract_primitives(buffer, 2632, 1570, 21, ncols);

                    simdfunc::contract_primitives(buffer, 2674, 498, 21, ncols);

                    simdfunc::contract_primitives(buffer, 2716, 1591, 21, ncols);

                    simdfunc::contract_primitives(buffer, 2758, 1612, 21, ncols);

                    simdfunc::contract_primitives(buffer, 2800, 1633, 21, ncols);

                    simdfunc::contract_primitives(buffer, 2842, 603, 21, ncols);

                    simdfunc::contract_primitives(buffer, 2884, 1654, 28, ncols);

                    simdfunc::contract_primitives(buffer, 2940, 1682, 28, ncols);

                    simdfunc::contract_primitives(buffer, 2996, 1710, 28, ncols);

                    simdfunc::contract_primitives(buffer, 3052, 708, 28, ncols);

                    simdfunc::contract_primitives(buffer, 3108, 1738, 28, ncols);

                    simdfunc::contract_primitives(buffer, 3164, 1766, 28, ncols);

                    simdfunc::contract_primitives(buffer, 3220, 1794, 28, ncols);

                    simdfunc::contract_primitives(buffer, 3276, 820, 28, ncols);

                    simdfunc::contract_primitives(buffer, 3332, 1822, 36, ncols);

                    simdfunc::contract_primitives(buffer, 3404, 1858, 36, ncols);

                    simdfunc::contract_primitives(buffer, 3476, 1894, 36, ncols);

                    simdfunc::contract_primitives(buffer, 3548, 932, 36, ncols);

                    simdfunc::contract_primitives(buffer, 3620, 1930, 36, ncols);

                    simdfunc::contract_primitives(buffer, 3692, 1966, 36, ncols);

                    simdfunc::contract_primitives(buffer, 3764, 2002, 36, ncols);

                    simdfunc::contract_primitives(buffer, 3836, 1040, 36, ncols);

                    simdfunc::contract_primitives(buffer, 3908, 2038, 45, ncols);

                    simdfunc::contract_primitives(buffer, 3998, 2083, 45, ncols);

                    simdfunc::contract_primitives(buffer, 4088, 2128, 45, ncols);

                    simdfunc::contract_primitives(buffer, 4178, 2173, 45, ncols);

                    simdfunc::contract_primitives(buffer, 4268, 2218, 45, ncols);

                    simdfunc::contract_primitives(buffer, 4358, 2263, 45, ncols);
                }
            }
        }

        simdtrf::transform_s_inner(buffer, 2323, 2308, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2353, 2338, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2383, 2368, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2413, 2398, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2443, 2428, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2473, 2458, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2503, 2488, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2533, 2518, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2569, 2548, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2611, 2590, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2653, 2632, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2695, 2674, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2737, 2716, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2779, 2758, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2821, 2800, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2863, 2842, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2912, 2884, 28, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2968, 2940, 28, 1, nmax);

        simdtrf::transform_s_inner(buffer, 3024, 2996, 28, 1, nmax);

        simdtrf::transform_s_inner(buffer, 3080, 3052, 28, 1, nmax);

        simdtrf::transform_s_inner(buffer, 3136, 3108, 28, 1, nmax);

        simdtrf::transform_s_inner(buffer, 3192, 3164, 28, 1, nmax);

        simdtrf::transform_s_inner(buffer, 3248, 3220, 28, 1, nmax);

        simdtrf::transform_s_inner(buffer, 3304, 3276, 28, 1, nmax);

        simdtrf::transform_s_inner(buffer, 3368, 3332, 36, 1, nmax);

        simdtrf::transform_s_inner(buffer, 3440, 3404, 36, 1, nmax);

        simdtrf::transform_s_inner(buffer, 3512, 3476, 36, 1, nmax);

        simdtrf::transform_s_inner(buffer, 3584, 3548, 36, 1, nmax);

        simdtrf::transform_s_inner(buffer, 3656, 3620, 36, 1, nmax);

        simdtrf::transform_s_inner(buffer, 3728, 3692, 36, 1, nmax);

        simdtrf::transform_s_inner(buffer, 3800, 3764, 36, 1, nmax);

        simdtrf::transform_s_inner(buffer, 3872, 3836, 36, 1, nmax);

        simdtrf::transform_s_inner(buffer, 3953, 3908, 45, 1, nmax);

        simdtrf::transform_s_inner(buffer, 4043, 3998, 45, 1, nmax);

        simdtrf::transform_s_inner(buffer, 4133, 4088, 45, 1, nmax);

        simdtrf::transform_s_inner(buffer, 4223, 4178, 45, 1, nmax);

        simdtrf::transform_s_inner(buffer, 4313, 4268, 45, 1, nmax);

        simdtrf::transform_s_inner(buffer, 4403, 4358, 45, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 4448, 2323, 2413, 2569, 1, nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 4493, 2353, 2413, 2611, 1, nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 4538, 2383, 2413, 2653, 1, nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 4583, 2413, 2695, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 4628, 2443, 2533, 2737, 1, nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 4673, 2473, 2533, 2779, 1, nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 4718, 2503, 2533, 2821, 1, nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 4763, 2533, 2863, 1, nmax);

        simdtrf::compute_hrr_geom_010x_ph(buffer, coordinates, 4808, 2569, 2695, 2912, 1, nmax);

        simdtrf::compute_hrr_geom_010y_ph(buffer, coordinates, 4871, 2611, 2695, 2968, 1, nmax);

        simdtrf::compute_hrr_geom_010z_ph(buffer, coordinates, 4934, 2653, 2695, 3024, 1, nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 4997, 2695, 3080, 1, nmax);

        simdtrf::compute_hrr_geom_010x_ph(buffer, coordinates, 5060, 2737, 2863, 3136, 1, nmax);

        simdtrf::compute_hrr_geom_010y_ph(buffer, coordinates, 5123, 2779, 2863, 3192, 1, nmax);

        simdtrf::compute_hrr_geom_010z_ph(buffer, coordinates, 5186, 2821, 2863, 3248, 1, nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 5249, 2863, 3304, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pi(buffer, coordinates, 5312, 2912, 3080, 3368, 1, nmax);

        simdtrf::compute_hrr_geom_010y_pi(buffer, coordinates, 5396, 2968, 3080, 3440, 1, nmax);

        simdtrf::compute_hrr_geom_010z_pi(buffer, coordinates, 5480, 3024, 3080, 3512, 1, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 5564, 3080, 3584, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pi(buffer, coordinates, 5648, 3136, 3304, 3656, 1, nmax);

        simdtrf::compute_hrr_geom_010y_pi(buffer, coordinates, 5732, 3192, 3304, 3728, 1, nmax);

        simdtrf::compute_hrr_geom_010z_pi(buffer, coordinates, 5816, 3248, 3304, 3800, 1, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 5900, 3304, 3872, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pk(buffer, coordinates, 5984, 3368, 3584, 3953, 1, nmax);

        simdtrf::compute_hrr_geom_010y_pk(buffer, coordinates, 6092, 3440, 3584, 4043, 1, nmax);

        simdtrf::compute_hrr_geom_010z_pk(buffer, coordinates, 6200, 3512, 3584, 4133, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pk(buffer, coordinates, 6308, 3656, 3872, 4223, 1, nmax);

        simdtrf::compute_hrr_geom_010y_pk(buffer, coordinates, 6416, 3728, 3872, 4313, 1, nmax);

        simdtrf::compute_hrr_geom_010z_pk(buffer, coordinates, 6524, 3800, 3872, 4403, 1, nmax);

        simdtrf::compute_hrr_geom_010x_dg(buffer, coordinates, 6632, 4448, 4583, 4808, 1, nmax);

        simdtrf::compute_hrr_geom_010y_dg(buffer, coordinates, 6722, 4493, 4583, 4871, 1, nmax);

        simdtrf::compute_hrr_geom_010z_dg(buffer, coordinates, 6812, 4538, 4583, 4934, 1, nmax);

        simdtrf::compute_hrr_dg(buffer, coordinates, 6902, 4583, 4997, 1, nmax);

        simdtrf::compute_hrr_geom_010x_dg(buffer, coordinates, 6992, 4628, 4763, 5060, 1, nmax);

        simdtrf::compute_hrr_geom_010y_dg(buffer, coordinates, 7082, 4673, 4763, 5123, 1, nmax);

        simdtrf::compute_hrr_geom_010z_dg(buffer, coordinates, 7172, 4718, 4763, 5186, 1, nmax);

        simdtrf::compute_hrr_dg(buffer, coordinates, 7262, 4763, 5249, 1, nmax);

        simdtrf::compute_hrr_geom_010x_dh(buffer, coordinates, 7352, 4808, 4997, 5312, 1, nmax);

        simdtrf::compute_hrr_geom_010y_dh(buffer, coordinates, 7478, 4871, 4997, 5396, 1, nmax);

        simdtrf::compute_hrr_geom_010z_dh(buffer, coordinates, 7604, 4934, 4997, 5480, 1, nmax);

        simdtrf::compute_hrr_dh(buffer, coordinates, 7730, 4997, 5564, 1, nmax);

        simdtrf::compute_hrr_geom_010x_dh(buffer, coordinates, 7856, 5060, 5249, 5648, 1, nmax);

        simdtrf::compute_hrr_geom_010y_dh(buffer, coordinates, 7982, 5123, 5249, 5732, 1, nmax);

        simdtrf::compute_hrr_geom_010z_dh(buffer, coordinates, 8108, 5186, 5249, 5816, 1, nmax);

        simdtrf::compute_hrr_dh(buffer, coordinates, 8234, 5249, 5900, 1, nmax);

        simdtrf::compute_hrr_geom_010x_di(buffer, coordinates, 8360, 5312, 5564, 5984, 1, nmax);

        simdtrf::compute_hrr_geom_010y_di(buffer, coordinates, 8528, 5396, 5564, 6092, 1, nmax);

        simdtrf::compute_hrr_geom_010z_di(buffer, coordinates, 8696, 5480, 5564, 6200, 1, nmax);

        simdtrf::compute_hrr_geom_010x_di(buffer, coordinates, 8864, 5648, 5900, 6308, 1, nmax);

        simdtrf::compute_hrr_geom_010y_di(buffer, coordinates, 9032, 5732, 5900, 6416, 1, nmax);

        simdtrf::compute_hrr_geom_010z_di(buffer, coordinates, 9200, 5816, 5900, 6524, 1, nmax);

        simdtrf::compute_hrr_geom_010x_fg(buffer, coordinates, 9368, 6632, 6902, 7352, 1, nmax);

        simdtrf::compute_hrr_geom_010y_fg(buffer, coordinates, 9518, 6722, 6902, 7478, 1, nmax);

        simdtrf::compute_hrr_geom_010z_fg(buffer, coordinates, 9668, 6812, 6902, 7604, 1, nmax);

        simdtrf::compute_hrr_fg(buffer, coordinates, 9818, 6902, 7730, 1, nmax);

        simdtrf::compute_hrr_geom_010x_fg(buffer, coordinates, 9968, 6992, 7262, 7856, 1, nmax);

        simdtrf::compute_hrr_geom_010y_fg(buffer, coordinates, 10118, 7082, 7262, 7982, 1,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_fg(buffer, coordinates, 10268, 7172, 7262, 8108, 1,
                                          nmax);

        simdtrf::compute_hrr_fg(buffer, coordinates, 10418, 7262, 8234, 1, nmax);

        simdtrf::compute_hrr_geom_010x_fh(buffer, coordinates, 10568, 7352, 7730, 8360, 1,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_fh(buffer, coordinates, 10778, 7478, 7730, 8528, 1,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_fh(buffer, coordinates, 10988, 7604, 7730, 8696, 1,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_fh(buffer, coordinates, 11198, 7856, 8234, 8864, 1,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_fh(buffer, coordinates, 11408, 7982, 8234, 9032, 1,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_fh(buffer, coordinates, 11618, 8108, 8234, 9200, 1,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_gg(buffer, coordinates, 11828, 9368, 9818, 10568, 1,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_gg(buffer, coordinates, 12053, 9518, 9818, 10778, 1,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_gg(buffer, coordinates, 12278, 9668, 9818, 10988, 1,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_gg(buffer, coordinates, 12503, 9968, 10418, 11198, 1,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_gg(buffer, coordinates, 12728, 10118, 10418, 11408, 1,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_gg(buffer, coordinates, 12953, 10268, 10418, 11618, 1,
                                          nmax);

        simdtrf::transform_g_inner(buffer, 13178, 12503, 15, 1, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 13178, 9, nmax);

        simdtrf::transform_g_inner(buffer, 13178, 12728, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 81 * nvalues + n * npairs, nvalues, buffer, 13178, 9,
                                   nmax);

        simdtrf::transform_g_inner(buffer, 13178, 12953, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 162 * nvalues + n * npairs, nvalues, buffer, 13178,
                                   9, nmax);

        simdtrf::transform_g_inner(buffer, 13178, 11828, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 243 * nvalues + n * npairs, nvalues, buffer, 13178,
                                   9, nmax);

        simdtrf::transform_g_inner(buffer, 13178, 12053, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 324 * nvalues + n * npairs, nvalues, buffer, 13178,
                                   9, nmax);

        simdtrf::transform_g_inner(buffer, 13178, 12278, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 405 * nvalues + n * npairs, nvalues, buffer, 13178,
                                   9, nmax);
    }

    for (size_t m = 0; m < 486; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
