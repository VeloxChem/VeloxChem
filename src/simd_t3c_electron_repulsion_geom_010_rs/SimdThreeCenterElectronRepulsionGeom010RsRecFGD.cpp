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


#include "SimdThreeCenterElectronRepulsionGeom010RsRecFGD.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSLD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferDG.hpp"
#include "SimdTransferGeom010XDG.hpp"
#include "SimdTransferGeom010XDH.hpp"
#include "SimdTransferGeom010XFG.hpp"
#include "SimdTransferGeom010XPG.hpp"
#include "SimdTransferGeom010XPH.hpp"
#include "SimdTransferGeom010XPI.hpp"
#include "SimdTransferGeom010YDG.hpp"
#include "SimdTransferGeom010YDH.hpp"
#include "SimdTransferGeom010YFG.hpp"
#include "SimdTransferGeom010YPG.hpp"
#include "SimdTransferGeom010YPH.hpp"
#include "SimdTransferGeom010YPI.hpp"
#include "SimdTransferGeom010ZDG.hpp"
#include "SimdTransferGeom010ZDH.hpp"
#include "SimdTransferGeom010ZFG.hpp"
#include "SimdTransferGeom010ZPG.hpp"
#include "SimdTransferGeom010ZPH.hpp"
#include "SimdTransferGeom010ZPI.hpp"
#include "SimdTransferPG.hpp"
#include "SimdTransferPH.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformG.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_010_fgd_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_010_fgd_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 41346, 0, 0, dimensions);

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
        simdfunc::prepare_buffer(buffer, 41346, 14168, 7828, dimensions);

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

                    simdfunc::compute_full_t3c_erf_boys_function(buffer, coordinates, 6, 3, 10,
                                                                 ncols, fj, i * nprim_b + j, fq,
                                                                 omega);

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 18, 3, 10,
                                                             ncols, fj, i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 30, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 33, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 36, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 39, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 42, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 45, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 48, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 51, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 54, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 57, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 60, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 63, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 66, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 69, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 72, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 75, 0, 3, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 78, 0, 3, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 81, 0, 3, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 84, 0, 3, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 87, 0, 3, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 90, 0, 3, 7, 8,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 96, 0, 3, 8, 9,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 102, 0, 3, 9, 10,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 108, 0, 3, 10, 11,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 114, 0, 3, 11, 12,
                                                                       42, 45, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 120, 0, 3, 12, 13,
                                                                       45, 48, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 126, 0, 3, 13, 14,
                                                                       48, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 132, 0, 3, 14, 15,
                                                                       51, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 138, 0, 3, 15, 16,
                                                                       54, 57, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 144, 0, 3, 19, 20,
                                                                       60, 63, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 150, 0, 3, 20, 21,
                                                                       63, 66, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 156, 0, 3, 21, 22,
                                                                       66, 69, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 162, 0, 3, 22, 23,
                                                                       69, 72, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 168, 0, 3, 23, 24,
                                                                       72, 75, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 174, 0, 3, 24, 25,
                                                                       75, 78, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 180, 0, 3, 25, 26,
                                                                       78, 81, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 186, 0, 3, 26, 27,
                                                                       81, 84, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 192, 0, 3, 27, 28,
                                                                       84, 87, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 198, 0, 3, 30, 33,
                                                                       90, 96, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 208, 0, 3, 33, 36,
                                                                       96, 102, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 218, 0, 3, 36, 39,
                                                                       102, 108, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 228, 0, 3, 39, 42,
                                                                       108, 114, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 238, 0, 3, 42, 45,
                                                                       114, 120, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 248, 0, 3, 45, 48,
                                                                       120, 126, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 258, 0, 3, 48, 51,
                                                                       126, 132, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 268, 0, 3, 51, 54,
                                                                       132, 138, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 278, 0, 3, 60, 63,
                                                                       144, 150, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 288, 0, 3, 63, 66,
                                                                       150, 156, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 298, 0, 3, 66, 69,
                                                                       156, 162, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 308, 0, 3, 69, 72,
                                                                       162, 168, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 318, 0, 3, 72, 75,
                                                                       168, 174, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 328, 0, 3, 75, 78,
                                                                       174, 180, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 338, 0, 3, 78, 81,
                                                                       180, 186, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 348, 0, 3, 81, 84,
                                                                       186, 192, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 358, 0, 3, 90, 96,
                                                                       198, 208, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 373, 0, 3, 96,
                                                                       102, 208, 218, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 388, 0, 3, 102,
                                                                       108, 218, 228, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 403, 0, 3, 108,
                                                                       114, 228, 238, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 418, 0, 3, 114,
                                                                       120, 238, 248, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 433, 0, 3, 120,
                                                                       126, 248, 258, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 448, 0, 3, 126,
                                                                       132, 258, 268, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 463, 0, 3, 144,
                                                                       150, 278, 288, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 478, 0, 3, 150,
                                                                       156, 288, 298, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 493, 0, 3, 156,
                                                                       162, 298, 308, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 508, 0, 3, 162,
                                                                       168, 308, 318, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 523, 0, 3, 168,
                                                                       174, 318, 328, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 538, 0, 3, 174,
                                                                       180, 328, 338, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 553, 0, 3, 180,
                                                                       186, 338, 348, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 568, 0, 3, 198,
                                                                       208, 358, 373, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 589, 0, 3, 208,
                                                                       218, 373, 388, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 610, 0, 3, 218,
                                                                       228, 388, 403, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 631, 0, 3, 228,
                                                                       238, 403, 418, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 652, 0, 3, 238,
                                                                       248, 418, 433, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 673, 0, 3, 248,
                                                                       258, 433, 448, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 694, 0, 3, 278,
                                                                       288, 463, 478, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 715, 0, 3, 288,
                                                                       298, 478, 493, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 736, 0, 3, 298,
                                                                       308, 493, 508, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 757, 0, 3, 308,
                                                                       318, 508, 523, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 778, 0, 3, 318,
                                                                       328, 523, 538, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 799, 0, 3, 328,
                                                                       338, 538, 553, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 820, 0, 3, 358,
                                                                       373, 568, 589, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 848, 0, 3, 373,
                                                                       388, 589, 610, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 876, 0, 3, 388,
                                                                       403, 610, 631, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 904, 0, 3, 403,
                                                                       418, 631, 652, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 932, 0, 3, 418,
                                                                       433, 652, 673, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 960, 0, 3, 463,
                                                                       478, 694, 715, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 988, 0, 3, 478,
                                                                       493, 715, 736, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1016, 0, 3, 493,
                                                                       508, 736, 757, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1044, 0, 3, 508,
                                                                       523, 757, 778, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1072, 0, 3, 523,
                                                                       538, 778, 799, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1100, 0, 3, 568,
                                                                       589, 820, 848, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1136, 0, 3, 589,
                                                                       610, 848, 876, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1172, 0, 3, 610,
                                                                       631, 876, 904, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1208, 0, 3, 631,
                                                                       652, 904, 932, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1244, 0, 3, 694,
                                                                       715, 960, 988, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1280, 0, 3, 715,
                                                                       736, 988, 1016, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1316, 0, 3, 736,
                                                                       757, 1016, 1044, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1352, 0, 3, 757,
                                                                       778, 1044, 1072, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1388, 0, 3, 820,
                                                                       848, 1100, 1136, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1433, 0, 3, 848,
                                                                       876, 1136, 1172, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1478, 0, 3, 876,
                                                                       904, 1172, 1208, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1523, 0, 3, 960,
                                                                       988, 1244, 1280, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1568, 0, 3, 988,
                                                                       1016, 1280, 1316, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1613, 0, 3, 1016,
                                                                       1044, 1316, 1352, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1658, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1661, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1664, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1667, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1670, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1673, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1676, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1679, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1682, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1685, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1688, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1691, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1694, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1697, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1700, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1703, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1706, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1709, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1712, 3, 9, 36,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1721, 3, 10, 39,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1730, 3, 11, 42,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1739, 3, 12, 45,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1748, 3, 13, 48,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1757, 3, 14, 51,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1766, 3, 15, 54,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1775, 3, 16, 57,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1784, 3, 21, 66,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1793, 3, 22, 69,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1802, 3, 23, 72,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1811, 3, 24, 75,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1820, 3, 25, 78,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1829, 3, 26, 81,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1838, 3, 27, 84,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1847, 3, 28, 87,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1856, 3, 36, 102,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1874, 3, 39, 108,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1892, 3, 42, 114,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1910, 3, 45, 120,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1928, 3, 48, 126,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1946, 3, 51, 132,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1964, 3, 54, 138,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1982, 3, 66, 156,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2000, 3, 69, 162,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2018, 3, 72, 168,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2036, 3, 75, 174,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2054, 3, 78, 180,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2072, 3, 81, 186,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2090, 3, 84, 192,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2108, 3, 102, 218,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2138, 3, 108, 228,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2168, 3, 114, 238,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2198, 3, 120, 248,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2228, 3, 126, 258,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2258, 3, 132, 268,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2288, 3, 156, 298,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2318, 3, 162, 308,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2348, 3, 168, 318,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2378, 3, 174, 328,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2408, 3, 180, 338,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2438, 3, 186, 348,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2468, 3, 218, 388,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2513, 3, 228, 403,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2558, 3, 238, 418,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2603, 3, 248, 433,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2648, 3, 258, 448,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2693, 3, 298, 493,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2738, 3, 308, 508,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2783, 3, 318, 523,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2828, 3, 328, 538,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2873, 3, 338, 553,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2918, 3, 388, 610,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2981, 3, 403, 631,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3044, 3, 418, 652,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3107, 3, 433, 673,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3170, 3, 493, 736,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3233, 3, 508, 757,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3296, 3, 523, 778,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3359, 3, 538, 799,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3422, 3, 610, 876,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3506, 3, 631, 904,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3590, 3, 652, 932,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3674, 3, 736,
                                                                       1016, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3758, 3, 757,
                                                                       1044, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3842, 3, 778,
                                                                       1072, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 3926, 3, 876,
                                                                       1172, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4034, 3, 904,
                                                                       1208, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4142, 3, 1016,
                                                                       1316, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4250, 3, 1044,
                                                                       1352, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 4358, 3, 1172,
                                                                       1478, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 4493, 3, 1316,
                                                                       1613, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4628, 3, 7, 8,
                                                                       1658, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4634, 3, 8, 9,
                                                                       1661, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4640, 3, 9, 10,
                                                                       1664, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4646, 3, 10, 11,
                                                                       1667, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4652, 3, 11, 12,
                                                                       1670, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4658, 3, 12, 13,
                                                                       1673, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4664, 3, 13, 14,
                                                                       1676, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4670, 3, 14, 15,
                                                                       1679, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4676, 3, 15, 16,
                                                                       1682, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4682, 3, 19, 20,
                                                                       1685, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4688, 3, 20, 21,
                                                                       1688, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4694, 3, 21, 22,
                                                                       1691, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4700, 3, 22, 23,
                                                                       1694, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4706, 3, 23, 24,
                                                                       1697, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4712, 3, 24, 25,
                                                                       1700, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4718, 3, 25, 26,
                                                                       1703, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4724, 3, 26, 27,
                                                                       1706, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4730, 3, 27, 28,
                                                                       1709, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4736, 0, 3, 4628,
                                                                       1658, 4634, 1712, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4754, 0, 3, 4634,
                                                                       1661, 4640, 1721, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4772, 0, 3, 4640,
                                                                       1664, 4646, 1730, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4790, 0, 3, 4646,
                                                                       1667, 4652, 1739, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4808, 0, 3, 4652,
                                                                       1670, 4658, 1748, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4826, 0, 3, 4658,
                                                                       1673, 4664, 1757, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4844, 0, 3, 4664,
                                                                       1676, 4670, 1766, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4862, 0, 3, 4670,
                                                                       1679, 4676, 1775, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4880, 0, 3, 4682,
                                                                       1685, 4688, 1784, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4898, 0, 3, 4688,
                                                                       1688, 4694, 1793, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4916, 0, 3, 4694,
                                                                       1691, 4700, 1802, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4934, 0, 3, 4700,
                                                                       1694, 4706, 1811, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4952, 0, 3, 4706,
                                                                       1697, 4712, 1820, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4970, 0, 3, 4712,
                                                                       1700, 4718, 1829, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4988, 0, 3, 4718,
                                                                       1703, 4724, 1838, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5006, 0, 3, 4724,
                                                                       1706, 4730, 1847, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5024, 0, 3, 4736,
                                                                       1712, 4754, 90, 96, 1856,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5060, 0, 3, 4754,
                                                                       1721, 4772, 96, 102, 1874,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5096, 0, 3, 4772,
                                                                       1730, 4790, 102, 108,
                                                                       1892, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5132, 0, 3, 4790,
                                                                       1739, 4808, 108, 114,
                                                                       1910, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5168, 0, 3, 4808,
                                                                       1748, 4826, 114, 120,
                                                                       1928, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5204, 0, 3, 4826,
                                                                       1757, 4844, 120, 126,
                                                                       1946, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5240, 0, 3, 4844,
                                                                       1766, 4862, 126, 132,
                                                                       1964, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5276, 0, 3, 4880,
                                                                       1784, 4898, 144, 150,
                                                                       1982, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5312, 0, 3, 4898,
                                                                       1793, 4916, 150, 156,
                                                                       2000, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5348, 0, 3, 4916,
                                                                       1802, 4934, 156, 162,
                                                                       2018, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5384, 0, 3, 4934,
                                                                       1811, 4952, 162, 168,
                                                                       2036, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5420, 0, 3, 4952,
                                                                       1820, 4970, 168, 174,
                                                                       2054, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5456, 0, 3, 4970,
                                                                       1829, 4988, 174, 180,
                                                                       2072, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5492, 0, 3, 4988,
                                                                       1838, 5006, 180, 186,
                                                                       2090, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5528, 0, 3, 5024,
                                                                       1856, 5060, 198, 208,
                                                                       2108, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5588, 0, 3, 5060,
                                                                       1874, 5096, 208, 218,
                                                                       2138, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5648, 0, 3, 5096,
                                                                       1892, 5132, 218, 228,
                                                                       2168, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5708, 0, 3, 5132,
                                                                       1910, 5168, 228, 238,
                                                                       2198, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5768, 0, 3, 5168,
                                                                       1928, 5204, 238, 248,
                                                                       2228, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5828, 0, 3, 5204,
                                                                       1946, 5240, 248, 258,
                                                                       2258, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5888, 0, 3, 5276,
                                                                       1982, 5312, 278, 288,
                                                                       2288, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5948, 0, 3, 5312,
                                                                       2000, 5348, 288, 298,
                                                                       2318, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6008, 0, 3, 5348,
                                                                       2018, 5384, 298, 308,
                                                                       2348, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6068, 0, 3, 5384,
                                                                       2036, 5420, 308, 318,
                                                                       2378, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6128, 0, 3, 5420,
                                                                       2054, 5456, 318, 328,
                                                                       2408, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6188, 0, 3, 5456,
                                                                       2072, 5492, 328, 338,
                                                                       2438, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6248, 0, 3, 5528,
                                                                       2108, 5588, 358, 373,
                                                                       2468, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6338, 0, 3, 5588,
                                                                       2138, 5648, 373, 388,
                                                                       2513, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6428, 0, 3, 5648,
                                                                       2168, 5708, 388, 403,
                                                                       2558, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6518, 0, 3, 5708,
                                                                       2198, 5768, 403, 418,
                                                                       2603, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6608, 0, 3, 5768,
                                                                       2228, 5828, 418, 433,
                                                                       2648, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6698, 0, 3, 5888,
                                                                       2288, 5948, 463, 478,
                                                                       2693, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6788, 0, 3, 5948,
                                                                       2318, 6008, 478, 493,
                                                                       2738, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6878, 0, 3, 6008,
                                                                       2348, 6068, 493, 508,
                                                                       2783, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6968, 0, 3, 6068,
                                                                       2378, 6128, 508, 523,
                                                                       2828, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7058, 0, 3, 6128,
                                                                       2408, 6188, 523, 538,
                                                                       2873, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 7148, 0, 3, 6248,
                                                                       2468, 6338, 568, 589,
                                                                       2918, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 7274, 0, 3, 6338,
                                                                       2513, 6428, 589, 610,
                                                                       2981, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 7400, 0, 3, 6428,
                                                                       2558, 6518, 610, 631,
                                                                       3044, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 7526, 0, 3, 6518,
                                                                       2603, 6608, 631, 652,
                                                                       3107, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 7652, 0, 3, 6698,
                                                                       2693, 6788, 694, 715,
                                                                       3170, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 7778, 0, 3, 6788,
                                                                       2738, 6878, 715, 736,
                                                                       3233, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 7904, 0, 3, 6878,
                                                                       2783, 6968, 736, 757,
                                                                       3296, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8030, 0, 3, 6968,
                                                                       2828, 7058, 757, 778,
                                                                       3359, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 8156, 0, 3, 7148,
                                                                       2918, 7274, 820, 848,
                                                                       3422, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 8324, 0, 3, 7274,
                                                                       2981, 7400, 848, 876,
                                                                       3506, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 8492, 0, 3, 7400,
                                                                       3044, 7526, 876, 904,
                                                                       3590, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 8660, 0, 3, 7652,
                                                                       3170, 7778, 960, 988,
                                                                       3674, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 8828, 0, 3, 7778,
                                                                       3233, 7904, 988, 1016,
                                                                       3758, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 8996, 0, 3, 7904,
                                                                       3296, 8030, 1016, 1044,
                                                                       3842, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 9164, 0, 3, 8156,
                                                                       3422, 8324, 1100, 1136,
                                                                       3926, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 9380, 0, 3, 8324,
                                                                       3506, 8492, 1136, 1172,
                                                                       4034, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 9596, 0, 3, 8660,
                                                                       3674, 8828, 1244, 1280,
                                                                       4142, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 9812, 0, 3, 8828,
                                                                       3758, 8996, 1280, 1316,
                                                                       4250, ncols, gamma, p,
                                                                       q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 10028, 0, 3, 9164,
                                                                       3926, 9380, 1388, 1433,
                                                                       4358, ncols, gamma, p,
                                                                       q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 10298, 0, 3, 9596,
                                                                       4142, 9812, 1523, 1568,
                                                                       4493, ncols, gamma, p,
                                                                       q);

                    simdgeo::geom_g_x(buffer, 10568, 5528, 7148, 1, 6, ncols, beta);

                    simdgeo::geom_g_y(buffer, 10658, 5528, 7148, 1, 6, ncols, beta);

                    simdgeo::geom_g_z(buffer, 10748, 5528, 7148, 1, 6, ncols, beta);

                    simdgeo::geom_g_x(buffer, 10838, 5888, 7652, 1, 6, ncols, beta);

                    simdgeo::geom_g_y(buffer, 10928, 5888, 7652, 1, 6, ncols, beta);

                    simdgeo::geom_g_z(buffer, 11018, 5888, 7652, 1, 6, ncols, beta);

                    simdgeo::geom_h_x(buffer, 11108, 6248, 8156, 1, 6, ncols, beta);

                    simdgeo::geom_h_y(buffer, 11234, 6248, 8156, 1, 6, ncols, beta);

                    simdgeo::geom_h_z(buffer, 11360, 6248, 8156, 1, 6, ncols, beta);

                    simdgeo::geom_h_x(buffer, 11486, 6698, 8660, 1, 6, ncols, beta);

                    simdgeo::geom_h_y(buffer, 11612, 6698, 8660, 1, 6, ncols, beta);

                    simdgeo::geom_h_z(buffer, 11738, 6698, 8660, 1, 6, ncols, beta);

                    simdgeo::geom_i_x(buffer, 11864, 7148, 9164, 1, 6, ncols, beta);

                    simdgeo::geom_i_y(buffer, 12032, 7148, 9164, 1, 6, ncols, beta);

                    simdgeo::geom_i_z(buffer, 12200, 7148, 9164, 1, 6, ncols, beta);

                    simdgeo::geom_i_x(buffer, 12368, 7652, 9596, 1, 6, ncols, beta);

                    simdgeo::geom_i_y(buffer, 12536, 7652, 9596, 1, 6, ncols, beta);

                    simdgeo::geom_i_z(buffer, 12704, 7652, 9596, 1, 6, ncols, beta);

                    simdgeo::geom_k_x(buffer, 12872, 8156, 10028, 1, 6, ncols, beta);

                    simdgeo::geom_k_y(buffer, 13088, 8156, 10028, 1, 6, ncols, beta);

                    simdgeo::geom_k_z(buffer, 13304, 8156, 10028, 1, 6, ncols, beta);

                    simdgeo::geom_k_x(buffer, 13520, 8660, 10298, 1, 6, ncols, beta);

                    simdgeo::geom_k_y(buffer, 13736, 8660, 10298, 1, 6, ncols, beta);

                    simdgeo::geom_k_z(buffer, 13952, 8660, 10298, 1, 6, ncols, beta);

                    simdfunc::contract_primitives(buffer, 14168, 10568, 90, ncols);

                    simdfunc::contract_primitives(buffer, 14333, 10658, 90, ncols);

                    simdfunc::contract_primitives(buffer, 14498, 10748, 90, ncols);

                    simdfunc::contract_primitives(buffer, 14663, 6248, 90, ncols);

                    simdfunc::contract_primitives(buffer, 14828, 10838, 90, ncols);

                    simdfunc::contract_primitives(buffer, 14993, 10928, 90, ncols);

                    simdfunc::contract_primitives(buffer, 15158, 11018, 90, ncols);

                    simdfunc::contract_primitives(buffer, 15323, 6698, 90, ncols);

                    simdfunc::contract_primitives(buffer, 15488, 11108, 126, ncols);

                    simdfunc::contract_primitives(buffer, 15719, 11234, 126, ncols);

                    simdfunc::contract_primitives(buffer, 15950, 11360, 126, ncols);

                    simdfunc::contract_primitives(buffer, 16181, 7148, 126, ncols);

                    simdfunc::contract_primitives(buffer, 16412, 11486, 126, ncols);

                    simdfunc::contract_primitives(buffer, 16643, 11612, 126, ncols);

                    simdfunc::contract_primitives(buffer, 16874, 11738, 126, ncols);

                    simdfunc::contract_primitives(buffer, 17105, 7652, 126, ncols);

                    simdfunc::contract_primitives(buffer, 17336, 11864, 168, ncols);

                    simdfunc::contract_primitives(buffer, 17644, 12032, 168, ncols);

                    simdfunc::contract_primitives(buffer, 17952, 12200, 168, ncols);

                    simdfunc::contract_primitives(buffer, 18260, 8156, 168, ncols);

                    simdfunc::contract_primitives(buffer, 18568, 12368, 168, ncols);

                    simdfunc::contract_primitives(buffer, 18876, 12536, 168, ncols);

                    simdfunc::contract_primitives(buffer, 19184, 12704, 168, ncols);

                    simdfunc::contract_primitives(buffer, 19492, 8660, 168, ncols);

                    simdfunc::contract_primitives(buffer, 19800, 12872, 216, ncols);

                    simdfunc::contract_primitives(buffer, 20196, 13088, 216, ncols);

                    simdfunc::contract_primitives(buffer, 20592, 13304, 216, ncols);

                    simdfunc::contract_primitives(buffer, 20988, 13520, 216, ncols);

                    simdfunc::contract_primitives(buffer, 21384, 13736, 216, ncols);

                    simdfunc::contract_primitives(buffer, 21780, 13952, 216, ncols);
                }
            }
        }

        simdtrf::transform_d_inner(buffer, 14258, 14168, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 14423, 14333, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 14588, 14498, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 14753, 14663, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 14918, 14828, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 15083, 14993, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 15248, 15158, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 15413, 15323, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 15614, 15488, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 15845, 15719, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 16076, 15950, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 16307, 16181, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 16538, 16412, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 16769, 16643, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 17000, 16874, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 17231, 17105, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 17504, 17336, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 17812, 17644, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 18120, 17952, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 18428, 18260, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 18736, 18568, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 19044, 18876, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 19352, 19184, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 19660, 19492, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 20016, 19800, 36, 1, nmax);

        simdtrf::transform_d_inner(buffer, 20412, 20196, 36, 1, nmax);

        simdtrf::transform_d_inner(buffer, 20808, 20592, 36, 1, nmax);

        simdtrf::transform_d_inner(buffer, 21204, 20988, 36, 1, nmax);

        simdtrf::transform_d_inner(buffer, 21600, 21384, 36, 1, nmax);

        simdtrf::transform_d_inner(buffer, 21996, 21780, 36, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 22176, 14258, 14753, 15614, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 22401, 14423, 14753, 15845, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 22626, 14588, 14753, 16076, 5,
                                          nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 22851, 14753, 16307, 5, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 23076, 14918, 15413, 16538, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 23301, 15083, 15413, 16769, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 23526, 15248, 15413, 17000, 5,
                                          nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 23751, 15413, 17231, 5, nmax);

        simdtrf::compute_hrr_geom_010x_ph(buffer, coordinates, 23976, 15614, 16307, 17504, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_ph(buffer, coordinates, 24291, 15845, 16307, 17812, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_ph(buffer, coordinates, 24606, 16076, 16307, 18120, 5,
                                          nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 24921, 16307, 18428, 5, nmax);

        simdtrf::compute_hrr_geom_010x_ph(buffer, coordinates, 25236, 16538, 17231, 18736, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_ph(buffer, coordinates, 25551, 16769, 17231, 19044, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_ph(buffer, coordinates, 25866, 17000, 17231, 19352, 5,
                                          nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 26181, 17231, 19660, 5, nmax);

        simdtrf::compute_hrr_geom_010x_pi(buffer, coordinates, 26496, 17504, 18428, 20016, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pi(buffer, coordinates, 26916, 17812, 18428, 20412, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pi(buffer, coordinates, 27336, 18120, 18428, 20808, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_pi(buffer, coordinates, 27756, 18736, 19660, 21204, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pi(buffer, coordinates, 28176, 19044, 19660, 21600, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pi(buffer, coordinates, 28596, 19352, 19660, 21996, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_dg(buffer, coordinates, 29016, 22176, 22851, 23976, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dg(buffer, coordinates, 29466, 22401, 22851, 24291, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dg(buffer, coordinates, 29916, 22626, 22851, 24606, 5,
                                          nmax);

        simdtrf::compute_hrr_dg(buffer, coordinates, 30366, 22851, 24921, 5, nmax);

        simdtrf::compute_hrr_geom_010x_dg(buffer, coordinates, 30816, 23076, 23751, 25236, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dg(buffer, coordinates, 31266, 23301, 23751, 25551, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dg(buffer, coordinates, 31716, 23526, 23751, 25866, 5,
                                          nmax);

        simdtrf::compute_hrr_dg(buffer, coordinates, 32166, 23751, 26181, 5, nmax);

        simdtrf::compute_hrr_geom_010x_dh(buffer, coordinates, 32616, 23976, 24921, 26496, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dh(buffer, coordinates, 33246, 24291, 24921, 26916, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dh(buffer, coordinates, 33876, 24606, 24921, 27336, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_dh(buffer, coordinates, 34506, 25236, 26181, 27756, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dh(buffer, coordinates, 35136, 25551, 26181, 28176, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dh(buffer, coordinates, 35766, 25866, 26181, 28596, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_fg(buffer, coordinates, 36396, 29016, 30366, 32616, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_fg(buffer, coordinates, 37146, 29466, 30366, 33246, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_fg(buffer, coordinates, 37896, 29916, 30366, 33876, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_fg(buffer, coordinates, 38646, 30816, 32166, 34506, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_fg(buffer, coordinates, 39396, 31266, 32166, 35136, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_fg(buffer, coordinates, 40146, 31716, 32166, 35766, 5,
                                          nmax);

        simdtrf::transform_g_inner(buffer, 40896, 38646, 10, 5, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 40896, 45, nmax);

        simdtrf::transform_g_inner(buffer, 40896, 39396, 10, 5, nmax);

        simdtrf::transform_f_outer(values + 315 * nvalues + n * npairs, nvalues, buffer, 40896,
                                   45, nmax);

        simdtrf::transform_g_inner(buffer, 40896, 40146, 10, 5, nmax);

        simdtrf::transform_f_outer(values + 630 * nvalues + n * npairs, nvalues, buffer, 40896,
                                   45, nmax);

        simdtrf::transform_g_inner(buffer, 40896, 36396, 10, 5, nmax);

        simdtrf::transform_f_outer(values + 945 * nvalues + n * npairs, nvalues, buffer, 40896,
                                   45, nmax);

        simdtrf::transform_g_inner(buffer, 40896, 37146, 10, 5, nmax);

        simdtrf::transform_f_outer(values + 1260 * nvalues + n * npairs, nvalues, buffer, 40896,
                                   45, nmax);

        simdtrf::transform_g_inner(buffer, 40896, 37896, 10, 5, nmax);

        simdtrf::transform_f_outer(values + 1575 * nvalues + n * npairs, nvalues, buffer, 40896,
                                   45, nmax);
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
