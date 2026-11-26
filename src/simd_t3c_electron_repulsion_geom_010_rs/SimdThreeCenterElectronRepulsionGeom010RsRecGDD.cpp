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


#include "SimdThreeCenterElectronRepulsionGeom010RsRecGDD.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSID.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
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
#include "SimdTransformG.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_010_gdd_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_010_gdd_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 37017, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1350 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 37017, 9968, 6284, dimensions);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1148, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1151, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1154, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1157, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1160, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1163, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1166, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1169, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1172, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1175, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1178, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1181, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1184, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1187, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1190, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1193, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1196, 3, 9, 34,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1205, 3, 10, 37,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1214, 3, 11, 40,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1223, 3, 12, 43,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1232, 3, 13, 46,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1241, 3, 14, 49,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1250, 3, 15, 52,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1259, 3, 20, 61,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1268, 3, 21, 64,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1277, 3, 22, 67,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1286, 3, 23, 70,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1295, 3, 24, 73,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1304, 3, 25, 76,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1313, 3, 26, 79,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1322, 3, 34, 94,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1340, 3, 37, 100,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1358, 3, 40, 106,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1376, 3, 43, 112,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1394, 3, 46, 118,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1412, 3, 49, 124,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1430, 3, 61, 142,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1448, 3, 64, 148,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1466, 3, 67, 154,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1484, 3, 70, 160,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1502, 3, 73, 166,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1520, 3, 76, 172,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1538, 3, 94, 198,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1568, 3, 100, 208,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1598, 3, 106, 218,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1628, 3, 112, 228,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1658, 3, 118, 238,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1688, 3, 142, 268,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1718, 3, 148, 278,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1748, 3, 154, 288,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1778, 3, 160, 298,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1808, 3, 166, 308,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1838, 3, 198, 348,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1883, 3, 208, 363,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1928, 3, 218, 378,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1973, 3, 228, 393,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2018, 3, 268, 438,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2063, 3, 278, 453,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2108, 3, 288, 468,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2153, 3, 298, 483,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2198, 3, 348, 540,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2261, 3, 363, 561,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2324, 3, 378, 582,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2387, 3, 438, 645,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2450, 3, 453, 666,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2513, 3, 468, 687,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2576, 3, 540, 764,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2660, 3, 561, 792,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2744, 3, 645, 876,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2828, 3, 666, 904,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 2912, 3, 764,
                                                                       1004, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 3020, 3, 876,
                                                                       1112, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3128, 3, 7, 8,
                                                                       1148, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3134, 3, 8, 9,
                                                                       1151, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3140, 3, 9, 10,
                                                                       1154, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3146, 3, 10, 11,
                                                                       1157, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3152, 3, 11, 12,
                                                                       1160, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3158, 3, 12, 13,
                                                                       1163, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3164, 3, 13, 14,
                                                                       1166, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3170, 3, 14, 15,
                                                                       1169, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3176, 3, 18, 19,
                                                                       1172, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3182, 3, 19, 20,
                                                                       1175, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3188, 3, 20, 21,
                                                                       1178, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3194, 3, 21, 22,
                                                                       1181, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3200, 3, 22, 23,
                                                                       1184, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3206, 3, 23, 24,
                                                                       1187, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3212, 3, 24, 25,
                                                                       1190, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3218, 3, 25, 26,
                                                                       1193, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3224, 0, 3, 3128,
                                                                       1148, 3134, 1196, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3242, 0, 3, 3134,
                                                                       1151, 3140, 1205, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3260, 0, 3, 3140,
                                                                       1154, 3146, 1214, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3278, 0, 3, 3146,
                                                                       1157, 3152, 1223, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3296, 0, 3, 3152,
                                                                       1160, 3158, 1232, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3314, 0, 3, 3158,
                                                                       1163, 3164, 1241, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3332, 0, 3, 3164,
                                                                       1166, 3170, 1250, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3350, 0, 3, 3176,
                                                                       1172, 3182, 1259, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3368, 0, 3, 3182,
                                                                       1175, 3188, 1268, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3386, 0, 3, 3188,
                                                                       1178, 3194, 1277, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3404, 0, 3, 3194,
                                                                       1181, 3200, 1286, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3422, 0, 3, 3200,
                                                                       1184, 3206, 1295, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3440, 0, 3, 3206,
                                                                       1187, 3212, 1304, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3458, 0, 3, 3212,
                                                                       1190, 3218, 1313, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3476, 0, 3, 3224,
                                                                       1196, 3242, 82, 88, 1322,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3512, 0, 3, 3242,
                                                                       1205, 3260, 88, 94, 1340,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3548, 0, 3, 3260,
                                                                       1214, 3278, 94, 100, 1358,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3584, 0, 3, 3278,
                                                                       1223, 3296, 100, 106,
                                                                       1376, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3620, 0, 3, 3296,
                                                                       1232, 3314, 106, 112,
                                                                       1394, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3656, 0, 3, 3314,
                                                                       1241, 3332, 112, 118,
                                                                       1412, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3692, 0, 3, 3350,
                                                                       1259, 3368, 130, 136,
                                                                       1430, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3728, 0, 3, 3368,
                                                                       1268, 3386, 136, 142,
                                                                       1448, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3764, 0, 3, 3386,
                                                                       1277, 3404, 142, 148,
                                                                       1466, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3800, 0, 3, 3404,
                                                                       1286, 3422, 148, 154,
                                                                       1484, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3836, 0, 3, 3422,
                                                                       1295, 3440, 154, 160,
                                                                       1502, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3872, 0, 3, 3440,
                                                                       1304, 3458, 160, 166,
                                                                       1520, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3908, 0, 3, 3476,
                                                                       1322, 3512, 178, 188,
                                                                       1538, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3968, 0, 3, 3512,
                                                                       1340, 3548, 188, 198,
                                                                       1568, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4028, 0, 3, 3548,
                                                                       1358, 3584, 198, 208,
                                                                       1598, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4088, 0, 3, 3584,
                                                                       1376, 3620, 208, 218,
                                                                       1628, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4148, 0, 3, 3620,
                                                                       1394, 3656, 218, 228,
                                                                       1658, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4208, 0, 3, 3692,
                                                                       1430, 3728, 248, 258,
                                                                       1688, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4268, 0, 3, 3728,
                                                                       1448, 3764, 258, 268,
                                                                       1718, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4328, 0, 3, 3764,
                                                                       1466, 3800, 268, 278,
                                                                       1748, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4388, 0, 3, 3800,
                                                                       1484, 3836, 278, 288,
                                                                       1778, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4448, 0, 3, 3836,
                                                                       1502, 3872, 288, 298,
                                                                       1808, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4508, 0, 3, 3908,
                                                                       1538, 3968, 318, 333,
                                                                       1838, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4598, 0, 3, 3968,
                                                                       1568, 4028, 333, 348,
                                                                       1883, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4688, 0, 3, 4028,
                                                                       1598, 4088, 348, 363,
                                                                       1928, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4778, 0, 3, 4088,
                                                                       1628, 4148, 363, 378,
                                                                       1973, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4868, 0, 3, 4208,
                                                                       1688, 4268, 408, 423,
                                                                       2018, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4958, 0, 3, 4268,
                                                                       1718, 4328, 423, 438,
                                                                       2063, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5048, 0, 3, 4328,
                                                                       1748, 4388, 438, 453,
                                                                       2108, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5138, 0, 3, 4388,
                                                                       1778, 4448, 453, 468,
                                                                       2153, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 5228, 0, 3, 4508,
                                                                       1838, 4598, 498, 519,
                                                                       2198, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 5354, 0, 3, 4598,
                                                                       1883, 4688, 519, 540,
                                                                       2261, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 5480, 0, 3, 4688,
                                                                       1928, 4778, 540, 561,
                                                                       2324, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 5606, 0, 3, 4868,
                                                                       2018, 4958, 603, 624,
                                                                       2387, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 5732, 0, 3, 4958,
                                                                       2063, 5048, 624, 645,
                                                                       2450, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 5858, 0, 3, 5048,
                                                                       2108, 5138, 645, 666,
                                                                       2513, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 5984, 0, 3, 5228,
                                                                       2198, 5354, 708, 736,
                                                                       2576, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 6152, 0, 3, 5354,
                                                                       2261, 5480, 736, 764,
                                                                       2660, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 6320, 0, 3, 5606,
                                                                       2387, 5732, 820, 848,
                                                                       2744, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 6488, 0, 3, 5732,
                                                                       2450, 5858, 848, 876,
                                                                       2828, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 6656, 0, 3, 5984,
                                                                       2576, 6152, 932, 968,
                                                                       2912, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 6872, 0, 3, 6320,
                                                                       2744, 6488, 1040, 1076,
                                                                       3020, ncols, gamma, p,
                                                                       q);

                    simdgeo::geom_d_x(buffer, 7088, 3224, 3908, 1, 6, ncols, beta);

                    simdgeo::geom_d_y(buffer, 7124, 3224, 3908, 1, 6, ncols, beta);

                    simdgeo::geom_d_z(buffer, 7160, 3224, 3908, 1, 6, ncols, beta);

                    simdgeo::geom_d_x(buffer, 7196, 3350, 4208, 1, 6, ncols, beta);

                    simdgeo::geom_d_y(buffer, 7232, 3350, 4208, 1, 6, ncols, beta);

                    simdgeo::geom_d_z(buffer, 7268, 3350, 4208, 1, 6, ncols, beta);

                    simdgeo::geom_f_x(buffer, 7304, 3476, 4508, 1, 6, ncols, beta);

                    simdgeo::geom_f_y(buffer, 7364, 3476, 4508, 1, 6, ncols, beta);

                    simdgeo::geom_f_z(buffer, 7424, 3476, 4508, 1, 6, ncols, beta);

                    simdgeo::geom_f_x(buffer, 7484, 3692, 4868, 1, 6, ncols, beta);

                    simdgeo::geom_f_y(buffer, 7544, 3692, 4868, 1, 6, ncols, beta);

                    simdgeo::geom_f_z(buffer, 7604, 3692, 4868, 1, 6, ncols, beta);

                    simdgeo::geom_g_x(buffer, 7664, 3908, 5228, 1, 6, ncols, beta);

                    simdgeo::geom_g_y(buffer, 7754, 3908, 5228, 1, 6, ncols, beta);

                    simdgeo::geom_g_z(buffer, 7844, 3908, 5228, 1, 6, ncols, beta);

                    simdgeo::geom_g_x(buffer, 7934, 4208, 5606, 1, 6, ncols, beta);

                    simdgeo::geom_g_y(buffer, 8024, 4208, 5606, 1, 6, ncols, beta);

                    simdgeo::geom_g_z(buffer, 8114, 4208, 5606, 1, 6, ncols, beta);

                    simdgeo::geom_h_x(buffer, 8204, 4508, 5984, 1, 6, ncols, beta);

                    simdgeo::geom_h_y(buffer, 8330, 4508, 5984, 1, 6, ncols, beta);

                    simdgeo::geom_h_z(buffer, 8456, 4508, 5984, 1, 6, ncols, beta);

                    simdgeo::geom_h_x(buffer, 8582, 4868, 6320, 1, 6, ncols, beta);

                    simdgeo::geom_h_y(buffer, 8708, 4868, 6320, 1, 6, ncols, beta);

                    simdgeo::geom_h_z(buffer, 8834, 4868, 6320, 1, 6, ncols, beta);

                    simdgeo::geom_i_x(buffer, 8960, 5228, 6656, 1, 6, ncols, beta);

                    simdgeo::geom_i_y(buffer, 9128, 5228, 6656, 1, 6, ncols, beta);

                    simdgeo::geom_i_z(buffer, 9296, 5228, 6656, 1, 6, ncols, beta);

                    simdgeo::geom_i_x(buffer, 9464, 5606, 6872, 1, 6, ncols, beta);

                    simdgeo::geom_i_y(buffer, 9632, 5606, 6872, 1, 6, ncols, beta);

                    simdgeo::geom_i_z(buffer, 9800, 5606, 6872, 1, 6, ncols, beta);

                    simdfunc::contract_primitives(buffer, 9968, 7088, 36, ncols);

                    simdfunc::contract_primitives(buffer, 10034, 7124, 36, ncols);

                    simdfunc::contract_primitives(buffer, 10100, 7160, 36, ncols);

                    simdfunc::contract_primitives(buffer, 10166, 3476, 36, ncols);

                    simdfunc::contract_primitives(buffer, 10232, 7196, 36, ncols);

                    simdfunc::contract_primitives(buffer, 10298, 7232, 36, ncols);

                    simdfunc::contract_primitives(buffer, 10364, 7268, 36, ncols);

                    simdfunc::contract_primitives(buffer, 10430, 3692, 36, ncols);

                    simdfunc::contract_primitives(buffer, 10496, 7304, 60, ncols);

                    simdfunc::contract_primitives(buffer, 10606, 7364, 60, ncols);

                    simdfunc::contract_primitives(buffer, 10716, 7424, 60, ncols);

                    simdfunc::contract_primitives(buffer, 10826, 3908, 60, ncols);

                    simdfunc::contract_primitives(buffer, 10936, 7484, 60, ncols);

                    simdfunc::contract_primitives(buffer, 11046, 7544, 60, ncols);

                    simdfunc::contract_primitives(buffer, 11156, 7604, 60, ncols);

                    simdfunc::contract_primitives(buffer, 11266, 4208, 60, ncols);

                    simdfunc::contract_primitives(buffer, 11376, 7664, 90, ncols);

                    simdfunc::contract_primitives(buffer, 11541, 7754, 90, ncols);

                    simdfunc::contract_primitives(buffer, 11706, 7844, 90, ncols);

                    simdfunc::contract_primitives(buffer, 11871, 4508, 90, ncols);

                    simdfunc::contract_primitives(buffer, 12036, 7934, 90, ncols);

                    simdfunc::contract_primitives(buffer, 12201, 8024, 90, ncols);

                    simdfunc::contract_primitives(buffer, 12366, 8114, 90, ncols);

                    simdfunc::contract_primitives(buffer, 12531, 4868, 90, ncols);

                    simdfunc::contract_primitives(buffer, 12696, 8204, 126, ncols);

                    simdfunc::contract_primitives(buffer, 12927, 8330, 126, ncols);

                    simdfunc::contract_primitives(buffer, 13158, 8456, 126, ncols);

                    simdfunc::contract_primitives(buffer, 13389, 5228, 126, ncols);

                    simdfunc::contract_primitives(buffer, 13620, 8582, 126, ncols);

                    simdfunc::contract_primitives(buffer, 13851, 8708, 126, ncols);

                    simdfunc::contract_primitives(buffer, 14082, 8834, 126, ncols);

                    simdfunc::contract_primitives(buffer, 14313, 5606, 126, ncols);

                    simdfunc::contract_primitives(buffer, 14544, 8960, 168, ncols);

                    simdfunc::contract_primitives(buffer, 14852, 9128, 168, ncols);

                    simdfunc::contract_primitives(buffer, 15160, 9296, 168, ncols);

                    simdfunc::contract_primitives(buffer, 15468, 9464, 168, ncols);

                    simdfunc::contract_primitives(buffer, 15776, 9632, 168, ncols);

                    simdfunc::contract_primitives(buffer, 16084, 9800, 168, ncols);
                }
            }
        }

        simdtrf::transform_d_inner(buffer, 10004, 9968, 6, 1, nmax);

        simdtrf::transform_d_inner(buffer, 10070, 10034, 6, 1, nmax);

        simdtrf::transform_d_inner(buffer, 10136, 10100, 6, 1, nmax);

        simdtrf::transform_d_inner(buffer, 10202, 10166, 6, 1, nmax);

        simdtrf::transform_d_inner(buffer, 10268, 10232, 6, 1, nmax);

        simdtrf::transform_d_inner(buffer, 10334, 10298, 6, 1, nmax);

        simdtrf::transform_d_inner(buffer, 10400, 10364, 6, 1, nmax);

        simdtrf::transform_d_inner(buffer, 10466, 10430, 6, 1, nmax);

        simdtrf::transform_d_inner(buffer, 10556, 10496, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 10666, 10606, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 10776, 10716, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 10886, 10826, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 10996, 10936, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 11106, 11046, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 11216, 11156, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 11326, 11266, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 11466, 11376, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 11631, 11541, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 11796, 11706, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 11961, 11871, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 12126, 12036, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 12291, 12201, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 12456, 12366, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 12621, 12531, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 12822, 12696, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 13053, 12927, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 13284, 13158, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 13515, 13389, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 13746, 13620, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 13977, 13851, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 14208, 14082, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 14439, 14313, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 14712, 14544, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 15020, 14852, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 15328, 15160, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 15636, 15468, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 15944, 15776, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 16252, 16084, 28, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pd(buffer, coordinates, 16392, 10004, 10202, 10556, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pd(buffer, coordinates, 16482, 10070, 10202, 10666, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pd(buffer, coordinates, 16572, 10136, 10202, 10776, 5,
                                          nmax);

        simdtrf::compute_hrr_pd(buffer, coordinates, 16662, 10202, 10886, 5, nmax);

        simdtrf::compute_hrr_geom_010x_pd(buffer, coordinates, 16752, 10268, 10466, 10996, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pd(buffer, coordinates, 16842, 10334, 10466, 11106, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pd(buffer, coordinates, 16932, 10400, 10466, 11216, 5,
                                          nmax);

        simdtrf::compute_hrr_pd(buffer, coordinates, 17022, 10466, 11326, 5, nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 17112, 10556, 10886, 11466, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 17262, 10666, 10886, 11631, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 17412, 10776, 10886, 11796, 5,
                                          nmax);

        simdtrf::compute_hrr_pf(buffer, coordinates, 17562, 10886, 11961, 5, nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 17712, 10996, 11326, 12126, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 17862, 11106, 11326, 12291, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 18012, 11216, 11326, 12456, 5,
                                          nmax);

        simdtrf::compute_hrr_pf(buffer, coordinates, 18162, 11326, 12621, 5, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 18312, 11466, 11961, 12822, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 18537, 11631, 11961, 13053, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 18762, 11796, 11961, 13284, 5,
                                          nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 18987, 11961, 13515, 5, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 19212, 12126, 12621, 13746, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 19437, 12291, 12621, 13977, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 19662, 12456, 12621, 14208, 5,
                                          nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 19887, 12621, 14439, 5, nmax);

        simdtrf::compute_hrr_geom_010x_ph(buffer, coordinates, 20112, 12822, 13515, 14712, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_ph(buffer, coordinates, 20427, 13053, 13515, 15020, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_ph(buffer, coordinates, 20742, 13284, 13515, 15328, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_ph(buffer, coordinates, 21057, 13746, 14439, 15636, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_ph(buffer, coordinates, 21372, 13977, 14439, 15944, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_ph(buffer, coordinates, 21687, 14208, 14439, 16252, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_dd(buffer, coordinates, 22002, 16392, 16662, 17112, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dd(buffer, coordinates, 22182, 16482, 16662, 17262, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dd(buffer, coordinates, 22362, 16572, 16662, 17412, 5,
                                          nmax);

        simdtrf::compute_hrr_dd(buffer, coordinates, 22542, 16662, 17562, 5, nmax);

        simdtrf::compute_hrr_geom_010x_dd(buffer, coordinates, 22722, 16752, 17022, 17712, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dd(buffer, coordinates, 22902, 16842, 17022, 17862, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dd(buffer, coordinates, 23082, 16932, 17022, 18012, 5,
                                          nmax);

        simdtrf::compute_hrr_dd(buffer, coordinates, 23262, 17022, 18162, 5, nmax);

        simdtrf::compute_hrr_geom_010x_df(buffer, coordinates, 23442, 17112, 17562, 18312, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_df(buffer, coordinates, 23742, 17262, 17562, 18537, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_df(buffer, coordinates, 24042, 17412, 17562, 18762, 5,
                                          nmax);

        simdtrf::compute_hrr_df(buffer, coordinates, 24342, 17562, 18987, 5, nmax);

        simdtrf::compute_hrr_geom_010x_df(buffer, coordinates, 24642, 17712, 18162, 19212, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_df(buffer, coordinates, 24942, 17862, 18162, 19437, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_df(buffer, coordinates, 25242, 18012, 18162, 19662, 5,
                                          nmax);

        simdtrf::compute_hrr_df(buffer, coordinates, 25542, 18162, 19887, 5, nmax);

        simdtrf::compute_hrr_geom_010x_dg(buffer, coordinates, 25842, 18312, 18987, 20112, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dg(buffer, coordinates, 26292, 18537, 18987, 20427, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dg(buffer, coordinates, 26742, 18762, 18987, 20742, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_dg(buffer, coordinates, 27192, 19212, 19887, 21057, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dg(buffer, coordinates, 27642, 19437, 19887, 21372, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dg(buffer, coordinates, 28092, 19662, 19887, 21687, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_fd_out_of_second(buffer, coordinates, 28542, 22002, 22542,
                                                        23442, 5, nmax);

        simdtrf::compute_hrr_geom_010y_fd_out_of_second(buffer, coordinates, 28842, 22182, 22542,
                                                        23742, 5, nmax);

        simdtrf::compute_hrr_geom_010z_fd_out_of_second(buffer, coordinates, 29142, 22362, 22542,
                                                        24042, 5, nmax);

        simdtrf::compute_hrr_fd_out_of_second(buffer, coordinates, 29442, 22542, 24342, 5,
                                              nmax);

        simdtrf::compute_hrr_geom_010x_fd_out_of_second(buffer, coordinates, 29742, 22722, 23262,
                                                        24642, 5, nmax);

        simdtrf::compute_hrr_geom_010y_fd_out_of_second(buffer, coordinates, 30042, 22902, 23262,
                                                        24942, 5, nmax);

        simdtrf::compute_hrr_geom_010z_fd_out_of_second(buffer, coordinates, 30342, 23082, 23262,
                                                        25242, 5, nmax);

        simdtrf::compute_hrr_fd_out_of_second(buffer, coordinates, 30642, 23262, 25542, 5,
                                              nmax);

        simdtrf::compute_hrr_geom_010x_ff(buffer, coordinates, 30942, 23442, 24342, 25842, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_ff(buffer, coordinates, 31442, 23742, 24342, 26292, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_ff(buffer, coordinates, 31942, 24042, 24342, 26742, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_ff(buffer, coordinates, 32442, 24642, 25542, 27192, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_ff(buffer, coordinates, 32942, 24942, 25542, 27642, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_ff(buffer, coordinates, 33442, 25242, 25542, 28092, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_gd_out_of_second(buffer, coordinates, 33942, 28542, 29442,
                                                        30942, 5, nmax);

        simdtrf::compute_hrr_geom_010y_gd_out_of_second(buffer, coordinates, 34392, 28842, 29442,
                                                        31442, 5, nmax);

        simdtrf::compute_hrr_geom_010z_gd_out_of_second(buffer, coordinates, 34842, 29142, 29442,
                                                        31942, 5, nmax);

        simdtrf::compute_hrr_geom_010x_gd_out_of_second(buffer, coordinates, 35292, 29742, 30642,
                                                        32442, 5, nmax);

        simdtrf::compute_hrr_geom_010y_gd_out_of_second(buffer, coordinates, 35742, 30042, 30642,
                                                        32942, 5, nmax);

        simdtrf::compute_hrr_geom_010z_gd_out_of_second(buffer, coordinates, 36192, 30342, 30642,
                                                        33442, 5, nmax);

        simdtrf::transform_d_inner(buffer, 36642, 35292, 15, 5, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 36642, 25, nmax);

        simdtrf::transform_d_inner(buffer, 36642, 35742, 15, 5, nmax);

        simdtrf::transform_g_outer(values + 225 * nvalues + n * npairs, nvalues, buffer, 36642,
                                   25, nmax);

        simdtrf::transform_d_inner(buffer, 36642, 36192, 15, 5, nmax);

        simdtrf::transform_g_outer(values + 450 * nvalues + n * npairs, nvalues, buffer, 36642,
                                   25, nmax);

        simdtrf::transform_d_inner(buffer, 36642, 33942, 15, 5, nmax);

        simdtrf::transform_g_outer(values + 675 * nvalues + n * npairs, nvalues, buffer, 36642,
                                   25, nmax);

        simdtrf::transform_d_inner(buffer, 36642, 34392, 15, 5, nmax);

        simdtrf::transform_g_outer(values + 900 * nvalues + n * npairs, nvalues, buffer, 36642,
                                   25, nmax);

        simdtrf::transform_d_inner(buffer, 36642, 34842, 15, 5, nmax);

        simdtrf::transform_g_outer(values + 1125 * nvalues + n * npairs, nvalues, buffer, 36642,
                                   25, nmax);
    }

    for (size_t m = 0; m < 1350; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
