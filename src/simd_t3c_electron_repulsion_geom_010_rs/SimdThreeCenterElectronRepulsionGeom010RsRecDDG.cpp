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


#include "SimdThreeCenterElectronRepulsionGeom010RsRecDDG.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSHD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferGeom010XDD.hpp"
#include "SimdTransferGeom010XPD.hpp"
#include "SimdTransferGeom010XPF.hpp"
#include "SimdTransferGeom010YDD.hpp"
#include "SimdTransferGeom010YPD.hpp"
#include "SimdTransferGeom010YPF.hpp"
#include "SimdTransferGeom010ZDD.hpp"
#include "SimdTransferGeom010ZPD.hpp"
#include "SimdTransferGeom010ZPF.hpp"
#include "SimdTransferPD.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformG.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_010_ddg_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_010_ddg_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 24444, 0, 0, dimensions);

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
        simdfunc::prepare_buffer(buffer, 24444, 14082, 5097, dimensions);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 708, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 711, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 714, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 717, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 720, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 723, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 726, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 729, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 732, 3, 20, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 735, 3, 21, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 738, 3, 22, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 741, 3, 23, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 744, 3, 24, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 747, 3, 25, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 750, 3, 26, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 753, 3, 27, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 756, 3, 9, 34,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 765, 3, 10, 37,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 774, 3, 11, 40,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 783, 3, 12, 43,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 792, 3, 13, 46,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 801, 3, 14, 49,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 810, 3, 15, 52,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 819, 3, 20, 61,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 828, 3, 21, 64,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 837, 3, 22, 67,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 846, 3, 23, 70,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 855, 3, 24, 73,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 864, 3, 25, 76,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 873, 3, 26, 79,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 882, 3, 34, 94,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 900, 3, 37, 100,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 918, 3, 40, 106,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 936, 3, 43, 112,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 954, 3, 46, 118,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 972, 3, 49, 124,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 990, 3, 61, 142,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1008, 3, 64, 148,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1026, 3, 67, 154,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1044, 3, 70, 160,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1062, 3, 73, 166,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1080, 3, 76, 172,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1098, 3, 94, 198,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1128, 3, 100, 208,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1158, 3, 106, 218,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1188, 3, 112, 228,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1218, 3, 118, 238,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1248, 3, 142, 268,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1278, 3, 148, 278,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1308, 3, 154, 288,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1338, 3, 160, 298,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1368, 3, 166, 308,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1398, 3, 198, 348,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1443, 3, 208, 363,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1488, 3, 218, 378,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1533, 3, 228, 393,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1578, 3, 268, 438,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1623, 3, 278, 453,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1668, 3, 288, 468,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1713, 3, 298, 483,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1758, 3, 348, 540,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1821, 3, 363, 561,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1884, 3, 378, 582,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1947, 3, 438, 645,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2010, 3, 453, 666,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2073, 3, 468, 687,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2136, 3, 7, 8,
                                                                       708, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2142, 3, 8, 9,
                                                                       711, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2148, 3, 9, 10,
                                                                       714, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2154, 3, 10, 11,
                                                                       717, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2160, 3, 11, 12,
                                                                       720, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2166, 3, 12, 13,
                                                                       723, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2172, 3, 13, 14,
                                                                       726, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2178, 3, 14, 15,
                                                                       729, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2184, 3, 18, 19,
                                                                       732, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2190, 3, 19, 20,
                                                                       735, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2196, 3, 20, 21,
                                                                       738, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2202, 3, 21, 22,
                                                                       741, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2208, 3, 22, 23,
                                                                       744, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2214, 3, 23, 24,
                                                                       747, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2220, 3, 24, 25,
                                                                       750, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2226, 3, 25, 26,
                                                                       753, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2232, 0, 3, 2136,
                                                                       708, 2142, 756, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2250, 0, 3, 2142,
                                                                       711, 2148, 765, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2268, 0, 3, 2148,
                                                                       714, 2154, 774, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2286, 0, 3, 2154,
                                                                       717, 2160, 783, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2304, 0, 3, 2160,
                                                                       720, 2166, 792, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2322, 0, 3, 2166,
                                                                       723, 2172, 801, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2340, 0, 3, 2172,
                                                                       726, 2178, 810, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2358, 0, 3, 2184,
                                                                       732, 2190, 819, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2376, 0, 3, 2190,
                                                                       735, 2196, 828, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2394, 0, 3, 2196,
                                                                       738, 2202, 837, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2412, 0, 3, 2202,
                                                                       741, 2208, 846, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2430, 0, 3, 2208,
                                                                       744, 2214, 855, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2448, 0, 3, 2214,
                                                                       747, 2220, 864, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2466, 0, 3, 2220,
                                                                       750, 2226, 873, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2484, 0, 3, 2232,
                                                                       756, 2250, 82, 88, 882,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2520, 0, 3, 2250,
                                                                       765, 2268, 88, 94, 900,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2556, 0, 3, 2268,
                                                                       774, 2286, 94, 100, 918,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2592, 0, 3, 2286,
                                                                       783, 2304, 100, 106, 936,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2628, 0, 3, 2304,
                                                                       792, 2322, 106, 112, 954,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2664, 0, 3, 2322,
                                                                       801, 2340, 112, 118, 972,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2700, 0, 3, 2358,
                                                                       819, 2376, 130, 136, 990,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2736, 0, 3, 2376,
                                                                       828, 2394, 136, 142, 1008,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2772, 0, 3, 2394,
                                                                       837, 2412, 142, 148, 1026,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2808, 0, 3, 2412,
                                                                       846, 2430, 148, 154, 1044,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2844, 0, 3, 2430,
                                                                       855, 2448, 154, 160, 1062,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2880, 0, 3, 2448,
                                                                       864, 2466, 160, 166, 1080,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2916, 0, 3, 2484,
                                                                       882, 2520, 178, 188, 1098,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2976, 0, 3, 2520,
                                                                       900, 2556, 188, 198, 1128,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3036, 0, 3, 2556,
                                                                       918, 2592, 198, 208, 1158,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3096, 0, 3, 2592,
                                                                       936, 2628, 208, 218, 1188,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3156, 0, 3, 2628,
                                                                       954, 2664, 218, 228, 1218,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3216, 0, 3, 2700,
                                                                       990, 2736, 248, 258, 1248,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3276, 0, 3, 2736,
                                                                       1008, 2772, 258, 268,
                                                                       1278, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3336, 0, 3, 2772,
                                                                       1026, 2808, 268, 278,
                                                                       1308, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3396, 0, 3, 2808,
                                                                       1044, 2844, 278, 288,
                                                                       1338, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3456, 0, 3, 2844,
                                                                       1062, 2880, 288, 298,
                                                                       1368, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3516, 0, 3, 2916,
                                                                       1098, 2976, 318, 333,
                                                                       1398, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3606, 0, 3, 2976,
                                                                       1128, 3036, 333, 348,
                                                                       1443, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3696, 0, 3, 3036,
                                                                       1158, 3096, 348, 363,
                                                                       1488, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3786, 0, 3, 3096,
                                                                       1188, 3156, 363, 378,
                                                                       1533, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3876, 0, 3, 3216,
                                                                       1248, 3276, 408, 423,
                                                                       1578, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3966, 0, 3, 3276,
                                                                       1278, 3336, 423, 438,
                                                                       1623, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4056, 0, 3, 3336,
                                                                       1308, 3396, 438, 453,
                                                                       1668, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4146, 0, 3, 3396,
                                                                       1338, 3456, 453, 468,
                                                                       1713, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 4236, 0, 3, 3516,
                                                                       1398, 3606, 498, 519,
                                                                       1758, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 4362, 0, 3, 3606,
                                                                       1443, 3696, 519, 540,
                                                                       1821, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 4488, 0, 3, 3696,
                                                                       1488, 3786, 540, 561,
                                                                       1884, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 4614, 0, 3, 3876,
                                                                       1578, 3966, 603, 624,
                                                                       1947, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 4740, 0, 3, 3966,
                                                                       1623, 4056, 624, 645,
                                                                       2010, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 4866, 0, 3, 4056,
                                                                       1668, 4146, 645, 666,
                                                                       2073, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4992, 3, 708, 711,
                                                                       2148, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5002, 3, 711, 714,
                                                                       2154, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5012, 3, 714, 717,
                                                                       2160, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5022, 3, 717, 720,
                                                                       2166, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5032, 3, 720, 723,
                                                                       2172, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5042, 3, 723, 726,
                                                                       2178, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5052, 3, 732, 735,
                                                                       2196, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5062, 3, 735, 738,
                                                                       2202, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5072, 3, 738, 741,
                                                                       2208, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5082, 3, 741, 744,
                                                                       2214, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5092, 3, 744, 747,
                                                                       2220, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5102, 3, 747, 750,
                                                                       2226, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5112, 0, 3, 4992,
                                                                       2148, 5002, 2268, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5142, 0, 3, 5002,
                                                                       2154, 5012, 2286, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5172, 0, 3, 5012,
                                                                       2160, 5022, 2304, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5202, 0, 3, 5022,
                                                                       2166, 5032, 2322, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5232, 0, 3, 5032,
                                                                       2172, 5042, 2340, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5262, 0, 3, 5052,
                                                                       2196, 5062, 2394, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5292, 0, 3, 5062,
                                                                       2202, 5072, 2412, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5322, 0, 3, 5072,
                                                                       2208, 5082, 2430, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5352, 0, 3, 5082,
                                                                       2214, 5092, 2448, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5382, 0, 3, 5092,
                                                                       2220, 5102, 2466, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 5412, 0, 3, 5112,
                                                                       2268, 5142, 882, 900,
                                                                       2556, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 5472, 0, 3, 5142,
                                                                       2286, 5172, 900, 918,
                                                                       2592, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 5532, 0, 3, 5172,
                                                                       2304, 5202, 918, 936,
                                                                       2628, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 5592, 0, 3, 5202,
                                                                       2322, 5232, 936, 954,
                                                                       2664, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 5652, 0, 3, 5262,
                                                                       2394, 5292, 990, 1008,
                                                                       2772, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 5712, 0, 3, 5292,
                                                                       2412, 5322, 1008, 1026,
                                                                       2808, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 5772, 0, 3, 5322,
                                                                       2430, 5352, 1026, 1044,
                                                                       2844, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 5832, 0, 3, 5352,
                                                                       2448, 5382, 1044, 1062,
                                                                       2880, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 5892, 0, 3, 5412,
                                                                       2556, 5472, 1098, 1128,
                                                                       3036, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 5992, 0, 3, 5472,
                                                                       2592, 5532, 1128, 1158,
                                                                       3096, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 6092, 0, 3, 5532,
                                                                       2628, 5592, 1158, 1188,
                                                                       3156, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 6192, 0, 3, 5652,
                                                                       2772, 5712, 1248, 1278,
                                                                       3336, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 6292, 0, 3, 5712,
                                                                       2808, 5772, 1278, 1308,
                                                                       3396, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 6392, 0, 3, 5772,
                                                                       2844, 5832, 1308, 1338,
                                                                       3456, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 6492, 0, 3, 5892,
                                                                       3036, 5992, 1398, 1443,
                                                                       3696, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 6642, 0, 3, 5992,
                                                                       3096, 6092, 1443, 1488,
                                                                       3786, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 6792, 0, 3, 6192,
                                                                       3336, 6292, 1578, 1623,
                                                                       4056, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 6942, 0, 3, 6292,
                                                                       3396, 6392, 1623, 1668,
                                                                       4146, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 7092, 0, 3, 6492,
                                                                       3696, 6642, 1758, 1821,
                                                                       4488, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 7302, 0, 3, 6792,
                                                                       4056, 6942, 1947, 2010,
                                                                       4866, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7512, 3, 2136,
                                                                       2142, 4992, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7527, 3, 2142,
                                                                       2148, 5002, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7542, 3, 2148,
                                                                       2154, 5012, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7557, 3, 2154,
                                                                       2160, 5022, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7572, 3, 2160,
                                                                       2166, 5032, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7587, 3, 2166,
                                                                       2172, 5042, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7602, 3, 2184,
                                                                       2190, 5052, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7617, 3, 2190,
                                                                       2196, 5062, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7632, 3, 2196,
                                                                       2202, 5072, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7647, 3, 2202,
                                                                       2208, 5082, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7662, 3, 2208,
                                                                       2214, 5092, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7677, 3, 2214,
                                                                       2220, 5102, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 7692, 0, 3, 7512,
                                                                       4992, 7527, 2232, 2250,
                                                                       5112, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 7737, 0, 3, 7527,
                                                                       5002, 7542, 2250, 2268,
                                                                       5142, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 7782, 0, 3, 7542,
                                                                       5012, 7557, 2268, 2286,
                                                                       5172, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 7827, 0, 3, 7557,
                                                                       5022, 7572, 2286, 2304,
                                                                       5202, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 7872, 0, 3, 7572,
                                                                       5032, 7587, 2304, 2322,
                                                                       5232, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 7917, 0, 3, 7602,
                                                                       5052, 7617, 2358, 2376,
                                                                       5262, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 7962, 0, 3, 7617,
                                                                       5062, 7632, 2376, 2394,
                                                                       5292, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8007, 0, 3, 7632,
                                                                       5072, 7647, 2394, 2412,
                                                                       5322, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8052, 0, 3, 7647,
                                                                       5082, 7662, 2412, 2430,
                                                                       5352, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8097, 0, 3, 7662,
                                                                       5092, 7677, 2430, 2448,
                                                                       5382, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 8142, 0, 3, 7692,
                                                                       5112, 7737, 2484, 2520,
                                                                       5412, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 8232, 0, 3, 7737,
                                                                       5142, 7782, 2520, 2556,
                                                                       5472, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 8322, 0, 3, 7782,
                                                                       5172, 7827, 2556, 2592,
                                                                       5532, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 8412, 0, 3, 7827,
                                                                       5202, 7872, 2592, 2628,
                                                                       5592, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 8502, 0, 3, 7917,
                                                                       5262, 7962, 2700, 2736,
                                                                       5652, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 8592, 0, 3, 7962,
                                                                       5292, 8007, 2736, 2772,
                                                                       5712, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 8682, 0, 3, 8007,
                                                                       5322, 8052, 2772, 2808,
                                                                       5772, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 8772, 0, 3, 8052,
                                                                       5352, 8097, 2808, 2844,
                                                                       5832, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 8862, 0, 3, 8142,
                                                                       5412, 8232, 2916, 2976,
                                                                       5892, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 9012, 0, 3, 8232,
                                                                       5472, 8322, 2976, 3036,
                                                                       5992, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 9162, 0, 3, 8322,
                                                                       5532, 8412, 3036, 3096,
                                                                       6092, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 9312, 0, 3, 8502,
                                                                       5652, 8592, 3216, 3276,
                                                                       6192, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 9462, 0, 3, 8592,
                                                                       5712, 8682, 3276, 3336,
                                                                       6292, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 9612, 0, 3, 8682,
                                                                       5772, 8772, 3336, 3396,
                                                                       6392, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 9762, 0, 3, 8862,
                                                                       5892, 9012, 3516, 3606,
                                                                       6492, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 9987, 0, 3, 9012,
                                                                       5992, 9162, 3606, 3696,
                                                                       6642, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 10212, 0, 3, 9312,
                                                                       6192, 9462, 3876, 3966,
                                                                       6792, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 10437, 0, 3, 9462,
                                                                       6292, 9612, 3966, 4056,
                                                                       6942, ncols, gamma, p,
                                                                       q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 10662, 0, 3, 9762,
                                                                       6492, 9987, 4236, 4362,
                                                                       7092, ncols, gamma, p,
                                                                       q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 10977, 0, 3,
                                                                       10212, 6792, 10437, 4614,
                                                                       4740, 7302, ncols, gamma,
                                                                       p, q);

                    simdgeo::geom_d_x(buffer, 11292, 7692, 8862, 1, 15, ncols, beta);

                    simdgeo::geom_d_y(buffer, 11382, 7692, 8862, 1, 15, ncols, beta);

                    simdgeo::geom_d_z(buffer, 11472, 7692, 8862, 1, 15, ncols, beta);

                    simdgeo::geom_d_x(buffer, 11562, 7917, 9312, 1, 15, ncols, beta);

                    simdgeo::geom_d_y(buffer, 11652, 7917, 9312, 1, 15, ncols, beta);

                    simdgeo::geom_d_z(buffer, 11742, 7917, 9312, 1, 15, ncols, beta);

                    simdgeo::geom_f_x(buffer, 11832, 8142, 9762, 1, 15, ncols, beta);

                    simdgeo::geom_f_y(buffer, 11982, 8142, 9762, 1, 15, ncols, beta);

                    simdgeo::geom_f_z(buffer, 12132, 8142, 9762, 1, 15, ncols, beta);

                    simdgeo::geom_f_x(buffer, 12282, 8502, 10212, 1, 15, ncols, beta);

                    simdgeo::geom_f_y(buffer, 12432, 8502, 10212, 1, 15, ncols, beta);

                    simdgeo::geom_f_z(buffer, 12582, 8502, 10212, 1, 15, ncols, beta);

                    simdgeo::geom_g_x(buffer, 12732, 8862, 10662, 1, 15, ncols, beta);

                    simdgeo::geom_g_y(buffer, 12957, 8862, 10662, 1, 15, ncols, beta);

                    simdgeo::geom_g_z(buffer, 13182, 8862, 10662, 1, 15, ncols, beta);

                    simdgeo::geom_g_x(buffer, 13407, 9312, 10977, 1, 15, ncols, beta);

                    simdgeo::geom_g_y(buffer, 13632, 9312, 10977, 1, 15, ncols, beta);

                    simdgeo::geom_g_z(buffer, 13857, 9312, 10977, 1, 15, ncols, beta);

                    simdfunc::contract_primitives(buffer, 14082, 11292, 90, ncols);

                    simdfunc::contract_primitives(buffer, 14226, 11382, 90, ncols);

                    simdfunc::contract_primitives(buffer, 14370, 11472, 90, ncols);

                    simdfunc::contract_primitives(buffer, 14514, 8142, 90, ncols);

                    simdfunc::contract_primitives(buffer, 14658, 11562, 90, ncols);

                    simdfunc::contract_primitives(buffer, 14802, 11652, 90, ncols);

                    simdfunc::contract_primitives(buffer, 14946, 11742, 90, ncols);

                    simdfunc::contract_primitives(buffer, 15090, 8502, 90, ncols);

                    simdfunc::contract_primitives(buffer, 15234, 11832, 150, ncols);

                    simdfunc::contract_primitives(buffer, 15474, 11982, 150, ncols);

                    simdfunc::contract_primitives(buffer, 15714, 12132, 150, ncols);

                    simdfunc::contract_primitives(buffer, 15954, 8862, 150, ncols);

                    simdfunc::contract_primitives(buffer, 16194, 12282, 150, ncols);

                    simdfunc::contract_primitives(buffer, 16434, 12432, 150, ncols);

                    simdfunc::contract_primitives(buffer, 16674, 12582, 150, ncols);

                    simdfunc::contract_primitives(buffer, 16914, 9312, 150, ncols);

                    simdfunc::contract_primitives(buffer, 17154, 12732, 225, ncols);

                    simdfunc::contract_primitives(buffer, 17514, 12957, 225, ncols);

                    simdfunc::contract_primitives(buffer, 17874, 13182, 225, ncols);

                    simdfunc::contract_primitives(buffer, 18234, 13407, 225, ncols);

                    simdfunc::contract_primitives(buffer, 18594, 13632, 225, ncols);

                    simdfunc::contract_primitives(buffer, 18954, 13857, 225, ncols);
                }
            }
        }

        simdtrf::transform_g_inner(buffer, 14172, 14082, 6, 1, nmax);

        simdtrf::transform_g_inner(buffer, 14316, 14226, 6, 1, nmax);

        simdtrf::transform_g_inner(buffer, 14460, 14370, 6, 1, nmax);

        simdtrf::transform_g_inner(buffer, 14604, 14514, 6, 1, nmax);

        simdtrf::transform_g_inner(buffer, 14748, 14658, 6, 1, nmax);

        simdtrf::transform_g_inner(buffer, 14892, 14802, 6, 1, nmax);

        simdtrf::transform_g_inner(buffer, 15036, 14946, 6, 1, nmax);

        simdtrf::transform_g_inner(buffer, 15180, 15090, 6, 1, nmax);

        simdtrf::transform_g_inner(buffer, 15384, 15234, 10, 1, nmax);

        simdtrf::transform_g_inner(buffer, 15624, 15474, 10, 1, nmax);

        simdtrf::transform_g_inner(buffer, 15864, 15714, 10, 1, nmax);

        simdtrf::transform_g_inner(buffer, 16104, 15954, 10, 1, nmax);

        simdtrf::transform_g_inner(buffer, 16344, 16194, 10, 1, nmax);

        simdtrf::transform_g_inner(buffer, 16584, 16434, 10, 1, nmax);

        simdtrf::transform_g_inner(buffer, 16824, 16674, 10, 1, nmax);

        simdtrf::transform_g_inner(buffer, 17064, 16914, 10, 1, nmax);

        simdtrf::transform_g_inner(buffer, 17379, 17154, 15, 1, nmax);

        simdtrf::transform_g_inner(buffer, 17739, 17514, 15, 1, nmax);

        simdtrf::transform_g_inner(buffer, 18099, 17874, 15, 1, nmax);

        simdtrf::transform_g_inner(buffer, 18459, 18234, 15, 1, nmax);

        simdtrf::transform_g_inner(buffer, 18819, 18594, 15, 1, nmax);

        simdtrf::transform_g_inner(buffer, 19179, 18954, 15, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pd(buffer, coordinates, 19314, 14172, 14604, 15384, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pd(buffer, coordinates, 19476, 14316, 14604, 15624, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pd(buffer, coordinates, 19638, 14460, 14604, 15864, 9,
                                          nmax);

        simdtrf::compute_hrr_pd(buffer, coordinates, 19800, 14604, 16104, 9, nmax);

        simdtrf::compute_hrr_geom_010x_pd(buffer, coordinates, 19962, 14748, 15180, 16344, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pd(buffer, coordinates, 20124, 14892, 15180, 16584, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pd(buffer, coordinates, 20286, 15036, 15180, 16824, 9,
                                          nmax);

        simdtrf::compute_hrr_pd(buffer, coordinates, 20448, 15180, 17064, 9, nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 20610, 15384, 16104, 17379, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 20880, 15624, 16104, 17739, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 21150, 15864, 16104, 18099, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 21420, 16344, 17064, 18459, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 21690, 16584, 17064, 18819, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 21960, 16824, 17064, 19179, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_dd(buffer, coordinates, 22230, 19314, 19800, 20610, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dd(buffer, coordinates, 22554, 19476, 19800, 20880, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dd(buffer, coordinates, 22878, 19638, 19800, 21150, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_dd(buffer, coordinates, 23202, 19962, 20448, 21420, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dd(buffer, coordinates, 23526, 20124, 20448, 21690, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dd(buffer, coordinates, 23850, 20286, 20448, 21960, 9,
                                          nmax);

        simdtrf::transform_d_inner(buffer, 24174, 23202, 6, 9, nmax);

        simdtrf::transform_d_outer(values + n * npairs, nvalues, buffer, 24174, 45, nmax);

        simdtrf::transform_d_inner(buffer, 24174, 23526, 6, 9, nmax);

        simdtrf::transform_d_outer(values + 225 * nvalues + n * npairs, nvalues, buffer, 24174,
                                   45, nmax);

        simdtrf::transform_d_inner(buffer, 24174, 23850, 6, 9, nmax);

        simdtrf::transform_d_outer(values + 450 * nvalues + n * npairs, nvalues, buffer, 24174,
                                   45, nmax);

        simdtrf::transform_d_inner(buffer, 24174, 22230, 6, 9, nmax);

        simdtrf::transform_d_outer(values + 675 * nvalues + n * npairs, nvalues, buffer, 24174,
                                   45, nmax);

        simdtrf::transform_d_inner(buffer, 24174, 22554, 6, 9, nmax);

        simdtrf::transform_d_outer(values + 900 * nvalues + n * npairs, nvalues, buffer, 24174,
                                   45, nmax);

        simdtrf::transform_d_inner(buffer, 24174, 22878, 6, 9, nmax);

        simdtrf::transform_d_outer(values + 1125 * nvalues + n * npairs, nvalues, buffer, 24174,
                                   45, nmax);
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
