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


#include "SimdThreeCenterElectronRepulsionGeom010RsRecGFF.hpp"

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
#include "SimdGeometryK1.hpp"
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
#include "SimdThreeCenterElectronRepulsionVrrRecSLD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferDF.hpp"
#include "SimdTransferDG.hpp"
#include "SimdTransferFF.hpp"
#include "SimdTransferGeom010XDF.hpp"
#include "SimdTransferGeom010XDG.hpp"
#include "SimdTransferGeom010XDH.hpp"
#include "SimdTransferGeom010XFF.hpp"
#include "SimdTransferGeom010XFG.hpp"
#include "SimdTransferGeom010XGF.hpp"
#include "SimdTransferGeom010XPF.hpp"
#include "SimdTransferGeom010XPG.hpp"
#include "SimdTransferGeom010XPH.hpp"
#include "SimdTransferGeom010XPI.hpp"
#include "SimdTransferGeom010YDF.hpp"
#include "SimdTransferGeom010YDG.hpp"
#include "SimdTransferGeom010YDH.hpp"
#include "SimdTransferGeom010YFF.hpp"
#include "SimdTransferGeom010YFG.hpp"
#include "SimdTransferGeom010YGF.hpp"
#include "SimdTransferGeom010YPF.hpp"
#include "SimdTransferGeom010YPG.hpp"
#include "SimdTransferGeom010YPH.hpp"
#include "SimdTransferGeom010YPI.hpp"
#include "SimdTransferGeom010ZDF.hpp"
#include "SimdTransferGeom010ZDG.hpp"
#include "SimdTransferGeom010ZDH.hpp"
#include "SimdTransferGeom010ZFF.hpp"
#include "SimdTransferGeom010ZFG.hpp"
#include "SimdTransferGeom010ZGF.hpp"
#include "SimdTransferGeom010ZPF.hpp"
#include "SimdTransferGeom010ZPG.hpp"
#include "SimdTransferGeom010ZPH.hpp"
#include "SimdTransferGeom010ZPI.hpp"
#include "SimdTransferPF.hpp"
#include "SimdTransferPG.hpp"
#include "SimdTransferPH.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformG.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_010_gff_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_010_gff_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 86631, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 2646 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 86631, 29012, 13484, dimensions);

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
                                                            4, 5, 6, 7, 8, 9, 10, 11}, ncols, fj,
                                                            i * nprim_b + j, fq, omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 18, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9, 10, 11}, ncols, fj,
                                                        i * nprim_b + j, fq);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1658, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1661, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1664, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1667, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1670, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1673, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1676, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1679, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1682, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1685, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1688, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1691, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1694, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1697, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1700, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1703, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1706, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1709, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1712, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1715, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1718, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1721, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1724, 3, 9, 36,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1733, 3, 10, 39,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1742, 3, 11, 42,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1751, 3, 12, 45,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1760, 3, 13, 48,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1769, 3, 14, 51,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1778, 3, 15, 54,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1787, 3, 16, 57,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1796, 3, 21, 66,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1805, 3, 22, 69,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1814, 3, 23, 72,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1823, 3, 24, 75,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1832, 3, 25, 78,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1841, 3, 26, 81,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1850, 3, 27, 84,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1859, 3, 28, 87,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1868, 3, 30, 90,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1886, 3, 33, 96,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1904, 3, 36, 102,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1922, 3, 39, 108,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1940, 3, 42, 114,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1958, 3, 45, 120,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1976, 3, 48, 126,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1994, 3, 51, 132,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2012, 3, 54, 138,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2030, 3, 60, 144,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2048, 3, 63, 150,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2066, 3, 66, 156,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2084, 3, 69, 162,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2102, 3, 72, 168,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2120, 3, 75, 174,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2138, 3, 78, 180,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2156, 3, 81, 186,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2174, 3, 84, 192,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2192, 3, 90, 198,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2222, 3, 96, 208,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2252, 3, 102, 218,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2282, 3, 108, 228,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2312, 3, 114, 238,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2342, 3, 120, 248,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2372, 3, 126, 258,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2402, 3, 132, 268,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2432, 3, 144, 278,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2462, 3, 150, 288,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2492, 3, 156, 298,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2522, 3, 162, 308,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2552, 3, 168, 318,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2582, 3, 174, 328,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2612, 3, 180, 338,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2642, 3, 186, 348,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2672, 3, 198, 358,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2717, 3, 208, 373,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2762, 3, 218, 388,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2807, 3, 228, 403,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2852, 3, 238, 418,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2897, 3, 248, 433,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2942, 3, 258, 448,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2987, 3, 278, 463,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3032, 3, 288, 478,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3077, 3, 298, 493,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3122, 3, 308, 508,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3167, 3, 318, 523,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3212, 3, 328, 538,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3257, 3, 338, 553,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3302, 3, 358, 568,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3365, 3, 373, 589,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3428, 3, 388, 610,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3491, 3, 403, 631,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3554, 3, 418, 652,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3617, 3, 433, 673,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3680, 3, 463, 694,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3743, 3, 478, 715,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3806, 3, 493, 736,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3869, 3, 508, 757,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3932, 3, 523, 778,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3995, 3, 538, 799,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4058, 3, 568, 820,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4142, 3, 589, 848,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4226, 3, 610, 876,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4310, 3, 631, 904,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4394, 3, 652, 932,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4478, 3, 694, 960,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4562, 3, 715, 988,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4646, 3, 736,
                                                                       1016, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4730, 3, 757,
                                                                       1044, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4814, 3, 778,
                                                                       1072, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4898, 3, 820,
                                                                       1100, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5006, 3, 848,
                                                                       1136, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5114, 3, 876,
                                                                       1172, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5222, 3, 904,
                                                                       1208, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5330, 3, 960,
                                                                       1244, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5438, 3, 988,
                                                                       1280, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5546, 3, 1016,
                                                                       1316, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5654, 3, 1044,
                                                                       1352, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5762, 3, 1100,
                                                                       1388, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5897, 3, 1136,
                                                                       1433, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6032, 3, 1172,
                                                                       1478, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6167, 3, 1244,
                                                                       1523, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6302, 3, 1280,
                                                                       1568, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6437, 3, 1316,
                                                                       1613, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6572, 3, 7, 8,
                                                                       1664, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6578, 3, 8, 9,
                                                                       1667, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6584, 3, 9, 10,
                                                                       1670, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6590, 3, 10, 11,
                                                                       1673, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6596, 3, 11, 12,
                                                                       1676, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6602, 3, 12, 13,
                                                                       1679, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6608, 3, 13, 14,
                                                                       1682, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6614, 3, 14, 15,
                                                                       1685, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6620, 3, 15, 16,
                                                                       1688, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6626, 3, 19, 20,
                                                                       1697, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6632, 3, 20, 21,
                                                                       1700, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6638, 3, 21, 22,
                                                                       1703, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6644, 3, 22, 23,
                                                                       1706, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6650, 3, 23, 24,
                                                                       1709, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6656, 3, 24, 25,
                                                                       1712, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6662, 3, 25, 26,
                                                                       1715, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6668, 3, 26, 27,
                                                                       1718, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6674, 3, 27, 28,
                                                                       1721, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6680, 0, 3, 6572,
                                                                       1664, 6578, 1724, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6698, 0, 3, 6578,
                                                                       1667, 6584, 1733, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6716, 0, 3, 6584,
                                                                       1670, 6590, 1742, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6734, 0, 3, 6590,
                                                                       1673, 6596, 1751, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6752, 0, 3, 6596,
                                                                       1676, 6602, 1760, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6770, 0, 3, 6602,
                                                                       1679, 6608, 1769, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6788, 0, 3, 6608,
                                                                       1682, 6614, 1778, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6806, 0, 3, 6614,
                                                                       1685, 6620, 1787, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6824, 0, 3, 6626,
                                                                       1697, 6632, 1796, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6842, 0, 3, 6632,
                                                                       1700, 6638, 1805, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6860, 0, 3, 6638,
                                                                       1703, 6644, 1814, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6878, 0, 3, 6644,
                                                                       1706, 6650, 1823, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6896, 0, 3, 6650,
                                                                       1709, 6656, 1832, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6914, 0, 3, 6656,
                                                                       1712, 6662, 1841, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6932, 0, 3, 6662,
                                                                       1715, 6668, 1850, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6950, 0, 3, 6668,
                                                                       1718, 6674, 1859, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6968, 0, 3, 6680,
                                                                       1724, 6698, 90, 96, 1904,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7004, 0, 3, 6698,
                                                                       1733, 6716, 96, 102, 1922,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7040, 0, 3, 6716,
                                                                       1742, 6734, 102, 108,
                                                                       1940, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7076, 0, 3, 6734,
                                                                       1751, 6752, 108, 114,
                                                                       1958, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7112, 0, 3, 6752,
                                                                       1760, 6770, 114, 120,
                                                                       1976, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7148, 0, 3, 6770,
                                                                       1769, 6788, 120, 126,
                                                                       1994, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7184, 0, 3, 6788,
                                                                       1778, 6806, 126, 132,
                                                                       2012, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7220, 0, 3, 6824,
                                                                       1796, 6842, 144, 150,
                                                                       2066, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7256, 0, 3, 6842,
                                                                       1805, 6860, 150, 156,
                                                                       2084, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7292, 0, 3, 6860,
                                                                       1814, 6878, 156, 162,
                                                                       2102, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7328, 0, 3, 6878,
                                                                       1823, 6896, 162, 168,
                                                                       2120, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7364, 0, 3, 6896,
                                                                       1832, 6914, 168, 174,
                                                                       2138, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7400, 0, 3, 6914,
                                                                       1841, 6932, 174, 180,
                                                                       2156, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7436, 0, 3, 6932,
                                                                       1850, 6950, 180, 186,
                                                                       2174, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7472, 0, 3, 6968,
                                                                       1904, 7004, 198, 208,
                                                                       2252, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7532, 0, 3, 7004,
                                                                       1922, 7040, 208, 218,
                                                                       2282, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7592, 0, 3, 7040,
                                                                       1940, 7076, 218, 228,
                                                                       2312, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7652, 0, 3, 7076,
                                                                       1958, 7112, 228, 238,
                                                                       2342, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7712, 0, 3, 7112,
                                                                       1976, 7148, 238, 248,
                                                                       2372, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7772, 0, 3, 7148,
                                                                       1994, 7184, 248, 258,
                                                                       2402, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7832, 0, 3, 7220,
                                                                       2066, 7256, 278, 288,
                                                                       2492, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7892, 0, 3, 7256,
                                                                       2084, 7292, 288, 298,
                                                                       2522, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7952, 0, 3, 7292,
                                                                       2102, 7328, 298, 308,
                                                                       2552, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8012, 0, 3, 7328,
                                                                       2120, 7364, 308, 318,
                                                                       2582, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8072, 0, 3, 7364,
                                                                       2138, 7400, 318, 328,
                                                                       2612, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8132, 0, 3, 7400,
                                                                       2156, 7436, 328, 338,
                                                                       2642, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8192, 0, 3, 7472,
                                                                       2252, 7532, 358, 373,
                                                                       2762, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8282, 0, 3, 7532,
                                                                       2282, 7592, 373, 388,
                                                                       2807, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8372, 0, 3, 7592,
                                                                       2312, 7652, 388, 403,
                                                                       2852, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8462, 0, 3, 7652,
                                                                       2342, 7712, 403, 418,
                                                                       2897, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8552, 0, 3, 7712,
                                                                       2372, 7772, 418, 433,
                                                                       2942, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8642, 0, 3, 7832,
                                                                       2492, 7892, 463, 478,
                                                                       3077, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8732, 0, 3, 7892,
                                                                       2522, 7952, 478, 493,
                                                                       3122, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8822, 0, 3, 7952,
                                                                       2552, 8012, 493, 508,
                                                                       3167, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8912, 0, 3, 8012,
                                                                       2582, 8072, 508, 523,
                                                                       3212, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9002, 0, 3, 8072,
                                                                       2612, 8132, 523, 538,
                                                                       3257, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9092, 0, 3, 8192,
                                                                       2762, 8282, 568, 589,
                                                                       3428, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9218, 0, 3, 8282,
                                                                       2807, 8372, 589, 610,
                                                                       3491, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9344, 0, 3, 8372,
                                                                       2852, 8462, 610, 631,
                                                                       3554, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9470, 0, 3, 8462,
                                                                       2897, 8552, 631, 652,
                                                                       3617, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9596, 0, 3, 8642,
                                                                       3077, 8732, 694, 715,
                                                                       3806, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9722, 0, 3, 8732,
                                                                       3122, 8822, 715, 736,
                                                                       3869, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9848, 0, 3, 8822,
                                                                       3167, 8912, 736, 757,
                                                                       3932, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9974, 0, 3, 8912,
                                                                       3212, 9002, 757, 778,
                                                                       3995, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 10100, 0, 3, 9092,
                                                                       3428, 9218, 820, 848,
                                                                       4226, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 10268, 0, 3, 9218,
                                                                       3491, 9344, 848, 876,
                                                                       4310, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 10436, 0, 3, 9344,
                                                                       3554, 9470, 876, 904,
                                                                       4394, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 10604, 0, 3, 9596,
                                                                       3806, 9722, 960, 988,
                                                                       4646, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 10772, 0, 3, 9722,
                                                                       3869, 9848, 988, 1016,
                                                                       4730, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 10940, 0, 3, 9848,
                                                                       3932, 9974, 1016, 1044,
                                                                       4814, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 11108, 0, 3,
                                                                       10100, 4226, 10268, 1100,
                                                                       1136, 5114, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 11324, 0, 3,
                                                                       10268, 4310, 10436, 1136,
                                                                       1172, 5222, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 11540, 0, 3,
                                                                       10604, 4646, 10772, 1244,
                                                                       1280, 5546, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 11756, 0, 3,
                                                                       10772, 4730, 10940, 1280,
                                                                       1316, 5654, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 11972, 0, 3,
                                                                       11108, 5114, 11324, 1388,
                                                                       1433, 6032, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 12242, 0, 3,
                                                                       11540, 5546, 11756, 1523,
                                                                       1568, 6437, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12512, 3, 1658,
                                                                       1661, 6572, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12522, 3, 1661,
                                                                       1664, 6578, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12532, 3, 1664,
                                                                       1667, 6584, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12542, 3, 1667,
                                                                       1670, 6590, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12552, 3, 1670,
                                                                       1673, 6596, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12562, 3, 1673,
                                                                       1676, 6602, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12572, 3, 1676,
                                                                       1679, 6608, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12582, 3, 1679,
                                                                       1682, 6614, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12592, 3, 1682,
                                                                       1685, 6620, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12602, 3, 1691,
                                                                       1694, 6626, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12612, 3, 1694,
                                                                       1697, 6632, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12622, 3, 1697,
                                                                       1700, 6638, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12632, 3, 1700,
                                                                       1703, 6644, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12642, 3, 1703,
                                                                       1706, 6650, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12652, 3, 1706,
                                                                       1709, 6656, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12662, 3, 1709,
                                                                       1712, 6662, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12672, 3, 1712,
                                                                       1715, 6668, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12682, 3, 1715,
                                                                       1718, 6674, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12692, 0, 3,
                                                                       12512, 6572, 12522, 6680,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12722, 0, 3,
                                                                       12522, 6578, 12532, 6698,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12752, 0, 3,
                                                                       12532, 6584, 12542, 6716,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12782, 0, 3,
                                                                       12542, 6590, 12552, 6734,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12812, 0, 3,
                                                                       12552, 6596, 12562, 6752,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12842, 0, 3,
                                                                       12562, 6602, 12572, 6770,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12872, 0, 3,
                                                                       12572, 6608, 12582, 6788,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12902, 0, 3,
                                                                       12582, 6614, 12592, 6806,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12932, 0, 3,
                                                                       12602, 6626, 12612, 6824,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12962, 0, 3,
                                                                       12612, 6632, 12622, 6842,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12992, 0, 3,
                                                                       12622, 6638, 12632, 6860,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13022, 0, 3,
                                                                       12632, 6644, 12642, 6878,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13052, 0, 3,
                                                                       12642, 6650, 12652, 6896,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13082, 0, 3,
                                                                       12652, 6656, 12662, 6914,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13112, 0, 3,
                                                                       12662, 6662, 12672, 6932,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13142, 0, 3,
                                                                       12672, 6668, 12682, 6950,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13172, 0, 3,
                                                                       12692, 6680, 12722, 1868,
                                                                       1886, 6968, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13232, 0, 3,
                                                                       12722, 6698, 12752, 1886,
                                                                       1904, 7004, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13292, 0, 3,
                                                                       12752, 6716, 12782, 1904,
                                                                       1922, 7040, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13352, 0, 3,
                                                                       12782, 6734, 12812, 1922,
                                                                       1940, 7076, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13412, 0, 3,
                                                                       12812, 6752, 12842, 1940,
                                                                       1958, 7112, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13472, 0, 3,
                                                                       12842, 6770, 12872, 1958,
                                                                       1976, 7148, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13532, 0, 3,
                                                                       12872, 6788, 12902, 1976,
                                                                       1994, 7184, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13592, 0, 3,
                                                                       12932, 6824, 12962, 2030,
                                                                       2048, 7220, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13652, 0, 3,
                                                                       12962, 6842, 12992, 2048,
                                                                       2066, 7256, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13712, 0, 3,
                                                                       12992, 6860, 13022, 2066,
                                                                       2084, 7292, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13772, 0, 3,
                                                                       13022, 6878, 13052, 2084,
                                                                       2102, 7328, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13832, 0, 3,
                                                                       13052, 6896, 13082, 2102,
                                                                       2120, 7364, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13892, 0, 3,
                                                                       13082, 6914, 13112, 2120,
                                                                       2138, 7400, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13952, 0, 3,
                                                                       13112, 6932, 13142, 2138,
                                                                       2156, 7436, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14012, 0, 3,
                                                                       13172, 6968, 13232, 2192,
                                                                       2222, 7472, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14112, 0, 3,
                                                                       13232, 7004, 13292, 2222,
                                                                       2252, 7532, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14212, 0, 3,
                                                                       13292, 7040, 13352, 2252,
                                                                       2282, 7592, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14312, 0, 3,
                                                                       13352, 7076, 13412, 2282,
                                                                       2312, 7652, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14412, 0, 3,
                                                                       13412, 7112, 13472, 2312,
                                                                       2342, 7712, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14512, 0, 3,
                                                                       13472, 7148, 13532, 2342,
                                                                       2372, 7772, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14612, 0, 3,
                                                                       13592, 7220, 13652, 2432,
                                                                       2462, 7832, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14712, 0, 3,
                                                                       13652, 7256, 13712, 2462,
                                                                       2492, 7892, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14812, 0, 3,
                                                                       13712, 7292, 13772, 2492,
                                                                       2522, 7952, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14912, 0, 3,
                                                                       13772, 7328, 13832, 2522,
                                                                       2552, 8012, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 15012, 0, 3,
                                                                       13832, 7364, 13892, 2552,
                                                                       2582, 8072, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 15112, 0, 3,
                                                                       13892, 7400, 13952, 2582,
                                                                       2612, 8132, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 15212, 0, 3,
                                                                       14012, 7472, 14112, 2672,
                                                                       2717, 8192, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 15362, 0, 3,
                                                                       14112, 7532, 14212, 2717,
                                                                       2762, 8282, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 15512, 0, 3,
                                                                       14212, 7592, 14312, 2762,
                                                                       2807, 8372, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 15662, 0, 3,
                                                                       14312, 7652, 14412, 2807,
                                                                       2852, 8462, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 15812, 0, 3,
                                                                       14412, 7712, 14512, 2852,
                                                                       2897, 8552, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 15962, 0, 3,
                                                                       14612, 7832, 14712, 2987,
                                                                       3032, 8642, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 16112, 0, 3,
                                                                       14712, 7892, 14812, 3032,
                                                                       3077, 8732, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 16262, 0, 3,
                                                                       14812, 7952, 14912, 3077,
                                                                       3122, 8822, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 16412, 0, 3,
                                                                       14912, 8012, 15012, 3122,
                                                                       3167, 8912, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 16562, 0, 3,
                                                                       15012, 8072, 15112, 3167,
                                                                       3212, 9002, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 16712, 0, 3,
                                                                       15212, 8192, 15362, 3302,
                                                                       3365, 9092, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 16922, 0, 3,
                                                                       15362, 8282, 15512, 3365,
                                                                       3428, 9218, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 17132, 0, 3,
                                                                       15512, 8372, 15662, 3428,
                                                                       3491, 9344, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 17342, 0, 3,
                                                                       15662, 8462, 15812, 3491,
                                                                       3554, 9470, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 17552, 0, 3,
                                                                       15962, 8642, 16112, 3680,
                                                                       3743, 9596, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 17762, 0, 3,
                                                                       16112, 8732, 16262, 3743,
                                                                       3806, 9722, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 17972, 0, 3,
                                                                       16262, 8822, 16412, 3806,
                                                                       3869, 9848, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 18182, 0, 3,
                                                                       16412, 8912, 16562, 3869,
                                                                       3932, 9974, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 18392, 0, 3,
                                                                       16712, 9092, 16922, 4058,
                                                                       4142, 10100, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 18672, 0, 3,
                                                                       16922, 9218, 17132, 4142,
                                                                       4226, 10268, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 18952, 0, 3,
                                                                       17132, 9344, 17342, 4226,
                                                                       4310, 10436, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 19232, 0, 3,
                                                                       17552, 9596, 17762, 4478,
                                                                       4562, 10604, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 19512, 0, 3,
                                                                       17762, 9722, 17972, 4562,
                                                                       4646, 10772, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 19792, 0, 3,
                                                                       17972, 9848, 18182, 4646,
                                                                       4730, 10940, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 20072, 0, 3,
                                                                       18392, 10100, 18672, 4898,
                                                                       5006, 11108, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 20432, 0, 3,
                                                                       18672, 10268, 18952, 5006,
                                                                       5114, 11324, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 20792, 0, 3,
                                                                       19232, 10604, 19512, 5330,
                                                                       5438, 11540, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 21152, 0, 3,
                                                                       19512, 10772, 19792, 5438,
                                                                       5546, 11756, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 21512, 0, 3,
                                                                       20072, 11108, 20432, 5762,
                                                                       5897, 11972, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 21962, 0, 3,
                                                                       20792, 11540, 21152, 6167,
                                                                       6302, 12242, ncols, gamma,
                                                                       p, q);

                    simdgeo::geom_f_x(buffer, 22412, 13172, 15212, 1, 10, ncols, beta);

                    simdgeo::geom_f_y(buffer, 22512, 13172, 15212, 1, 10, ncols, beta);

                    simdgeo::geom_f_z(buffer, 22612, 13172, 15212, 1, 10, ncols, beta);

                    simdgeo::geom_f_x(buffer, 22712, 13592, 15962, 1, 10, ncols, beta);

                    simdgeo::geom_f_y(buffer, 22812, 13592, 15962, 1, 10, ncols, beta);

                    simdgeo::geom_f_z(buffer, 22912, 13592, 15962, 1, 10, ncols, beta);

                    simdgeo::geom_g_x(buffer, 23012, 14012, 16712, 1, 10, ncols, beta);

                    simdgeo::geom_g_y(buffer, 23162, 14012, 16712, 1, 10, ncols, beta);

                    simdgeo::geom_g_z(buffer, 23312, 14012, 16712, 1, 10, ncols, beta);

                    simdgeo::geom_g_x(buffer, 23462, 14612, 17552, 1, 10, ncols, beta);

                    simdgeo::geom_g_y(buffer, 23612, 14612, 17552, 1, 10, ncols, beta);

                    simdgeo::geom_g_z(buffer, 23762, 14612, 17552, 1, 10, ncols, beta);

                    simdgeo::geom_h_x(buffer, 23912, 15212, 18392, 1, 10, ncols, beta);

                    simdgeo::geom_h_y(buffer, 24122, 15212, 18392, 1, 10, ncols, beta);

                    simdgeo::geom_h_z(buffer, 24332, 15212, 18392, 1, 10, ncols, beta);

                    simdgeo::geom_h_x(buffer, 24542, 15962, 19232, 1, 10, ncols, beta);

                    simdgeo::geom_h_y(buffer, 24752, 15962, 19232, 1, 10, ncols, beta);

                    simdgeo::geom_h_z(buffer, 24962, 15962, 19232, 1, 10, ncols, beta);

                    simdgeo::geom_i_x(buffer, 25172, 16712, 20072, 1, 10, ncols, beta);

                    simdgeo::geom_i_y(buffer, 25452, 16712, 20072, 1, 10, ncols, beta);

                    simdgeo::geom_i_z(buffer, 25732, 16712, 20072, 1, 10, ncols, beta);

                    simdgeo::geom_i_x(buffer, 26012, 17552, 20792, 1, 10, ncols, beta);

                    simdgeo::geom_i_y(buffer, 26292, 17552, 20792, 1, 10, ncols, beta);

                    simdgeo::geom_i_z(buffer, 26572, 17552, 20792, 1, 10, ncols, beta);

                    simdgeo::geom_k_x(buffer, 26852, 18392, 21512, 1, 10, ncols, beta);

                    simdgeo::geom_k_y(buffer, 27212, 18392, 21512, 1, 10, ncols, beta);

                    simdgeo::geom_k_z(buffer, 27572, 18392, 21512, 1, 10, ncols, beta);

                    simdgeo::geom_k_x(buffer, 27932, 19232, 21962, 1, 10, ncols, beta);

                    simdgeo::geom_k_y(buffer, 28292, 19232, 21962, 1, 10, ncols, beta);

                    simdgeo::geom_k_z(buffer, 28652, 19232, 21962, 1, 10, ncols, beta);

                    simdfunc::contract_primitives(buffer, 29012, 22412, 100, ncols);

                    simdfunc::contract_primitives(buffer, 29182, 22512, 100, ncols);

                    simdfunc::contract_primitives(buffer, 29352, 22612, 100, ncols);

                    simdfunc::contract_primitives(buffer, 29522, 14012, 100, ncols);

                    simdfunc::contract_primitives(buffer, 29692, 22712, 100, ncols);

                    simdfunc::contract_primitives(buffer, 29862, 22812, 100, ncols);

                    simdfunc::contract_primitives(buffer, 30032, 22912, 100, ncols);

                    simdfunc::contract_primitives(buffer, 30202, 14612, 100, ncols);

                    simdfunc::contract_primitives(buffer, 30372, 23012, 150, ncols);

                    simdfunc::contract_primitives(buffer, 30627, 23162, 150, ncols);

                    simdfunc::contract_primitives(buffer, 30882, 23312, 150, ncols);

                    simdfunc::contract_primitives(buffer, 31137, 15212, 150, ncols);

                    simdfunc::contract_primitives(buffer, 31392, 23462, 150, ncols);

                    simdfunc::contract_primitives(buffer, 31647, 23612, 150, ncols);

                    simdfunc::contract_primitives(buffer, 31902, 23762, 150, ncols);

                    simdfunc::contract_primitives(buffer, 32157, 15962, 150, ncols);

                    simdfunc::contract_primitives(buffer, 32412, 23912, 210, ncols);

                    simdfunc::contract_primitives(buffer, 32769, 24122, 210, ncols);

                    simdfunc::contract_primitives(buffer, 33126, 24332, 210, ncols);

                    simdfunc::contract_primitives(buffer, 33483, 16712, 210, ncols);

                    simdfunc::contract_primitives(buffer, 33840, 24542, 210, ncols);

                    simdfunc::contract_primitives(buffer, 34197, 24752, 210, ncols);

                    simdfunc::contract_primitives(buffer, 34554, 24962, 210, ncols);

                    simdfunc::contract_primitives(buffer, 34911, 17552, 210, ncols);

                    simdfunc::contract_primitives(buffer, 35268, 25172, 280, ncols);

                    simdfunc::contract_primitives(buffer, 35744, 25452, 280, ncols);

                    simdfunc::contract_primitives(buffer, 36220, 25732, 280, ncols);

                    simdfunc::contract_primitives(buffer, 36696, 18392, 280, ncols);

                    simdfunc::contract_primitives(buffer, 37172, 26012, 280, ncols);

                    simdfunc::contract_primitives(buffer, 37648, 26292, 280, ncols);

                    simdfunc::contract_primitives(buffer, 38124, 26572, 280, ncols);

                    simdfunc::contract_primitives(buffer, 38600, 19232, 280, ncols);

                    simdfunc::contract_primitives(buffer, 39076, 26852, 360, ncols);

                    simdfunc::contract_primitives(buffer, 39688, 27212, 360, ncols);

                    simdfunc::contract_primitives(buffer, 40300, 27572, 360, ncols);

                    simdfunc::contract_primitives(buffer, 40912, 27932, 360, ncols);

                    simdfunc::contract_primitives(buffer, 41524, 28292, 360, ncols);

                    simdfunc::contract_primitives(buffer, 42136, 28652, 360, ncols);
                }
            }
        }

        simdtrf::transform_f_inner(buffer, 29112, 29012, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 29282, 29182, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 29452, 29352, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 29622, 29522, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 29792, 29692, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 29962, 29862, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 30132, 30032, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 30302, 30202, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 30522, 30372, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 30777, 30627, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 31032, 30882, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 31287, 31137, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 31542, 31392, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 31797, 31647, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 32052, 31902, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 32307, 32157, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 32622, 32412, 21, 1, nmax);

        simdtrf::transform_f_inner(buffer, 32979, 32769, 21, 1, nmax);

        simdtrf::transform_f_inner(buffer, 33336, 33126, 21, 1, nmax);

        simdtrf::transform_f_inner(buffer, 33693, 33483, 21, 1, nmax);

        simdtrf::transform_f_inner(buffer, 34050, 33840, 21, 1, nmax);

        simdtrf::transform_f_inner(buffer, 34407, 34197, 21, 1, nmax);

        simdtrf::transform_f_inner(buffer, 34764, 34554, 21, 1, nmax);

        simdtrf::transform_f_inner(buffer, 35121, 34911, 21, 1, nmax);

        simdtrf::transform_f_inner(buffer, 35548, 35268, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 36024, 35744, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 36500, 36220, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 36976, 36696, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 37452, 37172, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 37928, 37648, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 38404, 38124, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 38880, 38600, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 39436, 39076, 36, 1, nmax);

        simdtrf::transform_f_inner(buffer, 40048, 39688, 36, 1, nmax);

        simdtrf::transform_f_inner(buffer, 40660, 40300, 36, 1, nmax);

        simdtrf::transform_f_inner(buffer, 41272, 40912, 36, 1, nmax);

        simdtrf::transform_f_inner(buffer, 41884, 41524, 36, 1, nmax);

        simdtrf::transform_f_inner(buffer, 42496, 42136, 36, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 42748, 29112, 29622, 30522, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 42958, 29282, 29622, 30777, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 43168, 29452, 29622, 31032, 7,
                                          nmax);

        simdtrf::compute_hrr_pf(buffer, coordinates, 43378, 29622, 31287, 7, nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 43588, 29792, 30302, 31542, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 43798, 29962, 30302, 31797, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 44008, 30132, 30302, 32052, 7,
                                          nmax);

        simdtrf::compute_hrr_pf(buffer, coordinates, 44218, 30302, 32307, 7, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 44428, 30522, 31287, 32622, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 44743, 30777, 31287, 32979, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 45058, 31032, 31287, 33336, 7,
                                          nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 45373, 31287, 33693, 7, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 45688, 31542, 32307, 34050, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 46003, 31797, 32307, 34407, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 46318, 32052, 32307, 34764, 7,
                                          nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 46633, 32307, 35121, 7, nmax);

        simdtrf::compute_hrr_geom_010x_ph(buffer, coordinates, 46948, 32622, 33693, 35548, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_ph(buffer, coordinates, 47389, 32979, 33693, 36024, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_ph(buffer, coordinates, 47830, 33336, 33693, 36500, 7,
                                          nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 48271, 33693, 36976, 7, nmax);

        simdtrf::compute_hrr_geom_010x_ph(buffer, coordinates, 48712, 34050, 35121, 37452, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_ph(buffer, coordinates, 49153, 34407, 35121, 37928, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_ph(buffer, coordinates, 49594, 34764, 35121, 38404, 7,
                                          nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 50035, 35121, 38880, 7, nmax);

        simdtrf::compute_hrr_geom_010x_pi(buffer, coordinates, 50476, 35548, 36976, 39436, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pi(buffer, coordinates, 51064, 36024, 36976, 40048, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pi(buffer, coordinates, 51652, 36500, 36976, 40660, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_pi(buffer, coordinates, 52240, 37452, 38880, 41272, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pi(buffer, coordinates, 52828, 37928, 38880, 41884, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pi(buffer, coordinates, 53416, 38404, 38880, 42496, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_df(buffer, coordinates, 54004, 42748, 43378, 44428, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_df(buffer, coordinates, 54424, 42958, 43378, 44743, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_df(buffer, coordinates, 54844, 43168, 43378, 45058, 7,
                                          nmax);

        simdtrf::compute_hrr_df(buffer, coordinates, 55264, 43378, 45373, 7, nmax);

        simdtrf::compute_hrr_geom_010x_df(buffer, coordinates, 55684, 43588, 44218, 45688, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_df(buffer, coordinates, 56104, 43798, 44218, 46003, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_df(buffer, coordinates, 56524, 44008, 44218, 46318, 7,
                                          nmax);

        simdtrf::compute_hrr_df(buffer, coordinates, 56944, 44218, 46633, 7, nmax);

        simdtrf::compute_hrr_geom_010x_dg(buffer, coordinates, 57364, 44428, 45373, 46948, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dg(buffer, coordinates, 57994, 44743, 45373, 47389, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dg(buffer, coordinates, 58624, 45058, 45373, 47830, 7,
                                          nmax);

        simdtrf::compute_hrr_dg(buffer, coordinates, 59254, 45373, 48271, 7, nmax);

        simdtrf::compute_hrr_geom_010x_dg(buffer, coordinates, 59884, 45688, 46633, 48712, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dg(buffer, coordinates, 60514, 46003, 46633, 49153, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dg(buffer, coordinates, 61144, 46318, 46633, 49594, 7,
                                          nmax);

        simdtrf::compute_hrr_dg(buffer, coordinates, 61774, 46633, 50035, 7, nmax);

        simdtrf::compute_hrr_geom_010x_dh(buffer, coordinates, 62404, 46948, 48271, 50476, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dh(buffer, coordinates, 63286, 47389, 48271, 51064, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dh(buffer, coordinates, 64168, 47830, 48271, 51652, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_dh(buffer, coordinates, 65050, 48712, 50035, 52240, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dh(buffer, coordinates, 65932, 49153, 50035, 52828, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dh(buffer, coordinates, 66814, 49594, 50035, 53416, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_ff(buffer, coordinates, 67696, 54004, 55264, 57364, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_ff(buffer, coordinates, 68396, 54424, 55264, 57994, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_ff(buffer, coordinates, 69096, 54844, 55264, 58624, 7,
                                          nmax);

        simdtrf::compute_hrr_ff(buffer, coordinates, 69796, 55264, 59254, 7, nmax);

        simdtrf::compute_hrr_geom_010x_ff(buffer, coordinates, 70496, 55684, 56944, 59884, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_ff(buffer, coordinates, 71196, 56104, 56944, 60514, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_ff(buffer, coordinates, 71896, 56524, 56944, 61144, 7,
                                          nmax);

        simdtrf::compute_hrr_ff(buffer, coordinates, 72596, 56944, 61774, 7, nmax);

        simdtrf::compute_hrr_geom_010x_fg(buffer, coordinates, 73296, 57364, 59254, 62404, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_fg(buffer, coordinates, 74346, 57994, 59254, 63286, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_fg(buffer, coordinates, 75396, 58624, 59254, 64168, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_fg(buffer, coordinates, 76446, 59884, 61774, 65050, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_fg(buffer, coordinates, 77496, 60514, 61774, 65932, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_fg(buffer, coordinates, 78546, 61144, 61774, 66814, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_gf_out_of_second(buffer, coordinates, 79596, 67696, 69796,
                                                        73296, 7, nmax);

        simdtrf::compute_hrr_geom_010y_gf_out_of_second(buffer, coordinates, 80646, 68396, 69796,
                                                        74346, 7, nmax);

        simdtrf::compute_hrr_geom_010z_gf_out_of_second(buffer, coordinates, 81696, 69096, 69796,
                                                        75396, 7, nmax);

        simdtrf::compute_hrr_geom_010x_gf_out_of_second(buffer, coordinates, 82746, 70496, 72596,
                                                        76446, 7, nmax);

        simdtrf::compute_hrr_geom_010y_gf_out_of_second(buffer, coordinates, 83796, 71196, 72596,
                                                        77496, 7, nmax);

        simdtrf::compute_hrr_geom_010z_gf_out_of_second(buffer, coordinates, 84846, 71896, 72596,
                                                        78546, 7, nmax);

        simdtrf::transform_f_inner(buffer, 85896, 82746, 15, 7, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 85896, 49, nmax);

        simdtrf::transform_f_inner(buffer, 85896, 83796, 15, 7, nmax);

        simdtrf::transform_g_outer(values + 441 * nvalues + n * npairs, nvalues, buffer, 85896,
                                   49, nmax);

        simdtrf::transform_f_inner(buffer, 85896, 84846, 15, 7, nmax);

        simdtrf::transform_g_outer(values + 882 * nvalues + n * npairs, nvalues, buffer, 85896,
                                   49, nmax);

        simdtrf::transform_f_inner(buffer, 85896, 79596, 15, 7, nmax);

        simdtrf::transform_g_outer(values + 1323 * nvalues + n * npairs, nvalues, buffer, 85896,
                                   49, nmax);

        simdtrf::transform_f_inner(buffer, 85896, 80646, 15, 7, nmax);

        simdtrf::transform_g_outer(values + 1764 * nvalues + n * npairs, nvalues, buffer, 85896,
                                   49, nmax);

        simdtrf::transform_f_inner(buffer, 85896, 81696, 15, 7, nmax);

        simdtrf::transform_g_outer(values + 2205 * nvalues + n * npairs, nvalues, buffer, 85896,
                                   49, nmax);
    }

    for (size_t m = 0; m < 2646; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
