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


#include "SimdThreeCenterElectronRepulsionGeom010RsRecGGD.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSMD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
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
#include "SimdTransformD.hpp"
#include "SimdTransformG.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_010_ggd_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_010_ggd_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 76503, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 2430 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 76503, 20408, 11545, dimensions);

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

                    simdfunc::compute_full_t3c_erf_boys_function(buffer, coordinates, 6, 3, 11,
                                                                 ncols, fj, i * nprim_b + j, fq,
                                                                 omega);

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 19, 3, 11,
                                                             ncols, fj, i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 32, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 35, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 38, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 41, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 44, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 47, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 50, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 53, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 56, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 59, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 62, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 65, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 68, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 71, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 74, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 77, 0, 3, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 80, 0, 3, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 83, 0, 3, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 86, 0, 3, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 89, 0, 3, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 92, 0, 3, 29, 30,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 95, 0, 3, 30, 31,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 98, 0, 3, 7, 8,
                                                                       32, 35, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 104, 0, 3, 8, 9,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 110, 0, 3, 9, 10,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 116, 0, 3, 10, 11,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 122, 0, 3, 11, 12,
                                                                       44, 47, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 128, 0, 3, 12, 13,
                                                                       47, 50, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 134, 0, 3, 13, 14,
                                                                       50, 53, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 140, 0, 3, 14, 15,
                                                                       53, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 146, 0, 3, 15, 16,
                                                                       56, 59, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 152, 0, 3, 16, 17,
                                                                       59, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 158, 0, 3, 20, 21,
                                                                       65, 68, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 164, 0, 3, 21, 22,
                                                                       68, 71, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 170, 0, 3, 22, 23,
                                                                       71, 74, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 176, 0, 3, 23, 24,
                                                                       74, 77, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 182, 0, 3, 24, 25,
                                                                       77, 80, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 188, 0, 3, 25, 26,
                                                                       80, 83, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 194, 0, 3, 26, 27,
                                                                       83, 86, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 200, 0, 3, 27, 28,
                                                                       86, 89, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 206, 0, 3, 28, 29,
                                                                       89, 92, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 212, 0, 3, 29, 30,
                                                                       92, 95, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 218, 0, 3, 32, 35,
                                                                       98, 104, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 228, 0, 3, 35, 38,
                                                                       104, 110, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 238, 0, 3, 38, 41,
                                                                       110, 116, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 248, 0, 3, 41, 44,
                                                                       116, 122, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 258, 0, 3, 44, 47,
                                                                       122, 128, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 268, 0, 3, 47, 50,
                                                                       128, 134, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 278, 0, 3, 50, 53,
                                                                       134, 140, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 288, 0, 3, 53, 56,
                                                                       140, 146, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 298, 0, 3, 56, 59,
                                                                       146, 152, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 308, 0, 3, 65, 68,
                                                                       158, 164, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 318, 0, 3, 68, 71,
                                                                       164, 170, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 328, 0, 3, 71, 74,
                                                                       170, 176, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 338, 0, 3, 74, 77,
                                                                       176, 182, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 348, 0, 3, 77, 80,
                                                                       182, 188, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 358, 0, 3, 80, 83,
                                                                       188, 194, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 368, 0, 3, 83, 86,
                                                                       194, 200, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 378, 0, 3, 86, 89,
                                                                       200, 206, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 388, 0, 3, 89, 92,
                                                                       206, 212, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 398, 0, 3, 98,
                                                                       104, 218, 228, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 413, 0, 3, 104,
                                                                       110, 228, 238, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 428, 0, 3, 110,
                                                                       116, 238, 248, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 443, 0, 3, 116,
                                                                       122, 248, 258, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 458, 0, 3, 122,
                                                                       128, 258, 268, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 473, 0, 3, 128,
                                                                       134, 268, 278, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 488, 0, 3, 134,
                                                                       140, 278, 288, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 503, 0, 3, 140,
                                                                       146, 288, 298, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 518, 0, 3, 158,
                                                                       164, 308, 318, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 533, 0, 3, 164,
                                                                       170, 318, 328, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 548, 0, 3, 170,
                                                                       176, 328, 338, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 563, 0, 3, 176,
                                                                       182, 338, 348, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 578, 0, 3, 182,
                                                                       188, 348, 358, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 593, 0, 3, 188,
                                                                       194, 358, 368, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 608, 0, 3, 194,
                                                                       200, 368, 378, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 623, 0, 3, 200,
                                                                       206, 378, 388, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 638, 0, 3, 218,
                                                                       228, 398, 413, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 659, 0, 3, 228,
                                                                       238, 413, 428, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 680, 0, 3, 238,
                                                                       248, 428, 443, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 701, 0, 3, 248,
                                                                       258, 443, 458, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 722, 0, 3, 258,
                                                                       268, 458, 473, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 743, 0, 3, 268,
                                                                       278, 473, 488, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 764, 0, 3, 278,
                                                                       288, 488, 503, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 785, 0, 3, 308,
                                                                       318, 518, 533, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 806, 0, 3, 318,
                                                                       328, 533, 548, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 827, 0, 3, 328,
                                                                       338, 548, 563, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 848, 0, 3, 338,
                                                                       348, 563, 578, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 869, 0, 3, 348,
                                                                       358, 578, 593, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 890, 0, 3, 358,
                                                                       368, 593, 608, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 911, 0, 3, 368,
                                                                       378, 608, 623, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 932, 0, 3, 398,
                                                                       413, 638, 659, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 960, 0, 3, 413,
                                                                       428, 659, 680, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 988, 0, 3, 428,
                                                                       443, 680, 701, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1016, 0, 3, 443,
                                                                       458, 701, 722, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1044, 0, 3, 458,
                                                                       473, 722, 743, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1072, 0, 3, 473,
                                                                       488, 743, 764, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1100, 0, 3, 518,
                                                                       533, 785, 806, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1128, 0, 3, 533,
                                                                       548, 806, 827, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1156, 0, 3, 548,
                                                                       563, 827, 848, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1184, 0, 3, 563,
                                                                       578, 848, 869, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1212, 0, 3, 578,
                                                                       593, 869, 890, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1240, 0, 3, 593,
                                                                       608, 890, 911, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1268, 0, 3, 638,
                                                                       659, 932, 960, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1304, 0, 3, 659,
                                                                       680, 960, 988, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1340, 0, 3, 680,
                                                                       701, 988, 1016, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1376, 0, 3, 701,
                                                                       722, 1016, 1044, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1412, 0, 3, 722,
                                                                       743, 1044, 1072, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1448, 0, 3, 785,
                                                                       806, 1100, 1128, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1484, 0, 3, 806,
                                                                       827, 1128, 1156, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1520, 0, 3, 827,
                                                                       848, 1156, 1184, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1556, 0, 3, 848,
                                                                       869, 1184, 1212, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1592, 0, 3, 869,
                                                                       890, 1212, 1240, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1628, 0, 3, 932,
                                                                       960, 1268, 1304, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1673, 0, 3, 960,
                                                                       988, 1304, 1340, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1718, 0, 3, 988,
                                                                       1016, 1340, 1376, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1763, 0, 3, 1016,
                                                                       1044, 1376, 1412, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1808, 0, 3, 1100,
                                                                       1128, 1448, 1484, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1853, 0, 3, 1128,
                                                                       1156, 1484, 1520, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1898, 0, 3, 1156,
                                                                       1184, 1520, 1556, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1943, 0, 3, 1184,
                                                                       1212, 1556, 1592, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1988, 0, 3, 1268,
                                                                       1304, 1628, 1673, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2043, 0, 3, 1304,
                                                                       1340, 1673, 1718, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2098, 0, 3, 1340,
                                                                       1376, 1718, 1763, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2153, 0, 3, 1448,
                                                                       1484, 1808, 1853, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2208, 0, 3, 1484,
                                                                       1520, 1853, 1898, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2263, 0, 3, 1520,
                                                                       1556, 1898, 1943, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2318, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2321, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2324, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2327, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2330, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2333, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2336, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2339, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2342, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2345, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2348, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2351, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2354, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2357, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2360, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2363, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2366, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2369, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2372, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2375, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2378, 3, 9, 38,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2387, 3, 10, 41,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2396, 3, 11, 44,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2405, 3, 12, 47,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2414, 3, 13, 50,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2423, 3, 14, 53,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2432, 3, 15, 56,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2441, 3, 16, 59,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2450, 3, 17, 62,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2459, 3, 22, 71,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2468, 3, 23, 74,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2477, 3, 24, 77,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2486, 3, 25, 80,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2495, 3, 26, 83,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2504, 3, 27, 86,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2513, 3, 28, 89,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2522, 3, 29, 92,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2531, 3, 30, 95,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2540, 3, 38, 110,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2558, 3, 41, 116,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2576, 3, 44, 122,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2594, 3, 47, 128,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2612, 3, 50, 134,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2630, 3, 53, 140,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2648, 3, 56, 146,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2666, 3, 59, 152,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2684, 3, 71, 170,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2702, 3, 74, 176,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2720, 3, 77, 182,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2738, 3, 80, 188,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2756, 3, 83, 194,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2774, 3, 86, 200,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2792, 3, 89, 206,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2810, 3, 92, 212,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2828, 3, 110, 238,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2858, 3, 116, 248,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2888, 3, 122, 258,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2918, 3, 128, 268,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2948, 3, 134, 278,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2978, 3, 140, 288,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3008, 3, 146, 298,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3038, 3, 170, 328,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3068, 3, 176, 338,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3098, 3, 182, 348,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3128, 3, 188, 358,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3158, 3, 194, 368,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3188, 3, 200, 378,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3218, 3, 206, 388,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3248, 3, 238, 428,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3293, 3, 248, 443,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3338, 3, 258, 458,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3383, 3, 268, 473,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3428, 3, 278, 488,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3473, 3, 288, 503,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3518, 3, 328, 548,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3563, 3, 338, 563,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3608, 3, 348, 578,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3653, 3, 358, 593,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3698, 3, 368, 608,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3743, 3, 378, 623,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3788, 3, 428, 680,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3851, 3, 443, 701,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3914, 3, 458, 722,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3977, 3, 473, 743,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4040, 3, 488, 764,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4103, 3, 548, 827,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4166, 3, 563, 848,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4229, 3, 578, 869,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4292, 3, 593, 890,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4355, 3, 608, 911,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4418, 3, 680, 988,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4502, 3, 701,
                                                                       1016, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4586, 3, 722,
                                                                       1044, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4670, 3, 743,
                                                                       1072, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4754, 3, 827,
                                                                       1156, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4838, 3, 848,
                                                                       1184, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4922, 3, 869,
                                                                       1212, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5006, 3, 890,
                                                                       1240, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5090, 3, 988,
                                                                       1340, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5198, 3, 1016,
                                                                       1376, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5306, 3, 1044,
                                                                       1412, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5414, 3, 1156,
                                                                       1520, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5522, 3, 1184,
                                                                       1556, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5630, 3, 1212,
                                                                       1592, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5738, 3, 1340,
                                                                       1718, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5873, 3, 1376,
                                                                       1763, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6008, 3, 1520,
                                                                       1898, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6143, 3, 1556,
                                                                       1943, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 6278, 3, 1718,
                                                                       2098, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 6443, 3, 1898,
                                                                       2263, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6608, 3, 7, 8,
                                                                       2318, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6614, 3, 8, 9,
                                                                       2321, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6620, 3, 9, 10,
                                                                       2324, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6626, 3, 10, 11,
                                                                       2327, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6632, 3, 11, 12,
                                                                       2330, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6638, 3, 12, 13,
                                                                       2333, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6644, 3, 13, 14,
                                                                       2336, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6650, 3, 14, 15,
                                                                       2339, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6656, 3, 15, 16,
                                                                       2342, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6662, 3, 16, 17,
                                                                       2345, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6668, 3, 20, 21,
                                                                       2348, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6674, 3, 21, 22,
                                                                       2351, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6680, 3, 22, 23,
                                                                       2354, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6686, 3, 23, 24,
                                                                       2357, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6692, 3, 24, 25,
                                                                       2360, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6698, 3, 25, 26,
                                                                       2363, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6704, 3, 26, 27,
                                                                       2366, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6710, 3, 27, 28,
                                                                       2369, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6716, 3, 28, 29,
                                                                       2372, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6722, 3, 29, 30,
                                                                       2375, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6728, 0, 3, 6608,
                                                                       2318, 6614, 2378, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6746, 0, 3, 6614,
                                                                       2321, 6620, 2387, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6764, 0, 3, 6620,
                                                                       2324, 6626, 2396, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6782, 0, 3, 6626,
                                                                       2327, 6632, 2405, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6800, 0, 3, 6632,
                                                                       2330, 6638, 2414, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6818, 0, 3, 6638,
                                                                       2333, 6644, 2423, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6836, 0, 3, 6644,
                                                                       2336, 6650, 2432, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6854, 0, 3, 6650,
                                                                       2339, 6656, 2441, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6872, 0, 3, 6656,
                                                                       2342, 6662, 2450, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6890, 0, 3, 6668,
                                                                       2348, 6674, 2459, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6908, 0, 3, 6674,
                                                                       2351, 6680, 2468, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6926, 0, 3, 6680,
                                                                       2354, 6686, 2477, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6944, 0, 3, 6686,
                                                                       2357, 6692, 2486, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6962, 0, 3, 6692,
                                                                       2360, 6698, 2495, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6980, 0, 3, 6698,
                                                                       2363, 6704, 2504, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6998, 0, 3, 6704,
                                                                       2366, 6710, 2513, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7016, 0, 3, 6710,
                                                                       2369, 6716, 2522, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7034, 0, 3, 6716,
                                                                       2372, 6722, 2531, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7052, 0, 3, 6728,
                                                                       2378, 6746, 98, 104, 2540,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7088, 0, 3, 6746,
                                                                       2387, 6764, 104, 110,
                                                                       2558, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7124, 0, 3, 6764,
                                                                       2396, 6782, 110, 116,
                                                                       2576, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7160, 0, 3, 6782,
                                                                       2405, 6800, 116, 122,
                                                                       2594, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7196, 0, 3, 6800,
                                                                       2414, 6818, 122, 128,
                                                                       2612, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7232, 0, 3, 6818,
                                                                       2423, 6836, 128, 134,
                                                                       2630, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7268, 0, 3, 6836,
                                                                       2432, 6854, 134, 140,
                                                                       2648, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7304, 0, 3, 6854,
                                                                       2441, 6872, 140, 146,
                                                                       2666, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7340, 0, 3, 6890,
                                                                       2459, 6908, 158, 164,
                                                                       2684, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7376, 0, 3, 6908,
                                                                       2468, 6926, 164, 170,
                                                                       2702, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7412, 0, 3, 6926,
                                                                       2477, 6944, 170, 176,
                                                                       2720, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7448, 0, 3, 6944,
                                                                       2486, 6962, 176, 182,
                                                                       2738, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7484, 0, 3, 6962,
                                                                       2495, 6980, 182, 188,
                                                                       2756, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7520, 0, 3, 6980,
                                                                       2504, 6998, 188, 194,
                                                                       2774, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7556, 0, 3, 6998,
                                                                       2513, 7016, 194, 200,
                                                                       2792, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7592, 0, 3, 7016,
                                                                       2522, 7034, 200, 206,
                                                                       2810, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7628, 0, 3, 7052,
                                                                       2540, 7088, 218, 228,
                                                                       2828, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7688, 0, 3, 7088,
                                                                       2558, 7124, 228, 238,
                                                                       2858, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7748, 0, 3, 7124,
                                                                       2576, 7160, 238, 248,
                                                                       2888, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7808, 0, 3, 7160,
                                                                       2594, 7196, 248, 258,
                                                                       2918, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7868, 0, 3, 7196,
                                                                       2612, 7232, 258, 268,
                                                                       2948, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7928, 0, 3, 7232,
                                                                       2630, 7268, 268, 278,
                                                                       2978, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7988, 0, 3, 7268,
                                                                       2648, 7304, 278, 288,
                                                                       3008, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8048, 0, 3, 7340,
                                                                       2684, 7376, 308, 318,
                                                                       3038, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8108, 0, 3, 7376,
                                                                       2702, 7412, 318, 328,
                                                                       3068, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8168, 0, 3, 7412,
                                                                       2720, 7448, 328, 338,
                                                                       3098, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8228, 0, 3, 7448,
                                                                       2738, 7484, 338, 348,
                                                                       3128, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8288, 0, 3, 7484,
                                                                       2756, 7520, 348, 358,
                                                                       3158, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8348, 0, 3, 7520,
                                                                       2774, 7556, 358, 368,
                                                                       3188, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8408, 0, 3, 7556,
                                                                       2792, 7592, 368, 378,
                                                                       3218, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8468, 0, 3, 7628,
                                                                       2828, 7688, 398, 413,
                                                                       3248, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8558, 0, 3, 7688,
                                                                       2858, 7748, 413, 428,
                                                                       3293, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8648, 0, 3, 7748,
                                                                       2888, 7808, 428, 443,
                                                                       3338, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8738, 0, 3, 7808,
                                                                       2918, 7868, 443, 458,
                                                                       3383, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8828, 0, 3, 7868,
                                                                       2948, 7928, 458, 473,
                                                                       3428, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8918, 0, 3, 7928,
                                                                       2978, 7988, 473, 488,
                                                                       3473, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9008, 0, 3, 8048,
                                                                       3038, 8108, 518, 533,
                                                                       3518, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9098, 0, 3, 8108,
                                                                       3068, 8168, 533, 548,
                                                                       3563, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9188, 0, 3, 8168,
                                                                       3098, 8228, 548, 563,
                                                                       3608, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9278, 0, 3, 8228,
                                                                       3128, 8288, 563, 578,
                                                                       3653, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9368, 0, 3, 8288,
                                                                       3158, 8348, 578, 593,
                                                                       3698, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9458, 0, 3, 8348,
                                                                       3188, 8408, 593, 608,
                                                                       3743, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9548, 0, 3, 8468,
                                                                       3248, 8558, 638, 659,
                                                                       3788, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9674, 0, 3, 8558,
                                                                       3293, 8648, 659, 680,
                                                                       3851, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9800, 0, 3, 8648,
                                                                       3338, 8738, 680, 701,
                                                                       3914, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9926, 0, 3, 8738,
                                                                       3383, 8828, 701, 722,
                                                                       3977, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10052, 0, 3, 8828,
                                                                       3428, 8918, 722, 743,
                                                                       4040, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10178, 0, 3, 9008,
                                                                       3518, 9098, 785, 806,
                                                                       4103, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10304, 0, 3, 9098,
                                                                       3563, 9188, 806, 827,
                                                                       4166, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10430, 0, 3, 9188,
                                                                       3608, 9278, 827, 848,
                                                                       4229, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10556, 0, 3, 9278,
                                                                       3653, 9368, 848, 869,
                                                                       4292, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10682, 0, 3, 9368,
                                                                       3698, 9458, 869, 890,
                                                                       4355, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 10808, 0, 3, 9548,
                                                                       3788, 9674, 932, 960,
                                                                       4418, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 10976, 0, 3, 9674,
                                                                       3851, 9800, 960, 988,
                                                                       4502, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11144, 0, 3, 9800,
                                                                       3914, 9926, 988, 1016,
                                                                       4586, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11312, 0, 3, 9926,
                                                                       3977, 10052, 1016, 1044,
                                                                       4670, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11480, 0, 3,
                                                                       10178, 4103, 10304, 1100,
                                                                       1128, 4754, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11648, 0, 3,
                                                                       10304, 4166, 10430, 1128,
                                                                       1156, 4838, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11816, 0, 3,
                                                                       10430, 4229, 10556, 1156,
                                                                       1184, 4922, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11984, 0, 3,
                                                                       10556, 4292, 10682, 1184,
                                                                       1212, 5006, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 12152, 0, 3,
                                                                       10808, 4418, 10976, 1268,
                                                                       1304, 5090, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 12368, 0, 3,
                                                                       10976, 4502, 11144, 1304,
                                                                       1340, 5198, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 12584, 0, 3,
                                                                       11144, 4586, 11312, 1340,
                                                                       1376, 5306, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 12800, 0, 3,
                                                                       11480, 4754, 11648, 1448,
                                                                       1484, 5414, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 13016, 0, 3,
                                                                       11648, 4838, 11816, 1484,
                                                                       1520, 5522, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 13232, 0, 3,
                                                                       11816, 4922, 11984, 1520,
                                                                       1556, 5630, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 13448, 0, 3,
                                                                       12152, 5090, 12368, 1628,
                                                                       1673, 5738, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 13718, 0, 3,
                                                                       12368, 5198, 12584, 1673,
                                                                       1718, 5873, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 13988, 0, 3,
                                                                       12800, 5414, 13016, 1808,
                                                                       1853, 6008, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 14258, 0, 3,
                                                                       13016, 5522, 13232, 1853,
                                                                       1898, 6143, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 14528, 0, 3,
                                                                       13448, 5738, 13718, 1988,
                                                                       2043, 6278, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 14858, 0, 3,
                                                                       13988, 6008, 14258, 2153,
                                                                       2208, 6443, ncols, gamma,
                                                                       p, q);

                    simdgeo::geom_g_x(buffer, 15188, 7628, 9548, 1, 6, ncols, beta);

                    simdgeo::geom_g_y(buffer, 15278, 7628, 9548, 1, 6, ncols, beta);

                    simdgeo::geom_g_z(buffer, 15368, 7628, 9548, 1, 6, ncols, beta);

                    simdgeo::geom_g_x(buffer, 15458, 8048, 10178, 1, 6, ncols, beta);

                    simdgeo::geom_g_y(buffer, 15548, 8048, 10178, 1, 6, ncols, beta);

                    simdgeo::geom_g_z(buffer, 15638, 8048, 10178, 1, 6, ncols, beta);

                    simdgeo::geom_h_x(buffer, 15728, 8468, 10808, 1, 6, ncols, beta);

                    simdgeo::geom_h_y(buffer, 15854, 8468, 10808, 1, 6, ncols, beta);

                    simdgeo::geom_h_z(buffer, 15980, 8468, 10808, 1, 6, ncols, beta);

                    simdgeo::geom_h_x(buffer, 16106, 9008, 11480, 1, 6, ncols, beta);

                    simdgeo::geom_h_y(buffer, 16232, 9008, 11480, 1, 6, ncols, beta);

                    simdgeo::geom_h_z(buffer, 16358, 9008, 11480, 1, 6, ncols, beta);

                    simdgeo::geom_i_x(buffer, 16484, 9548, 12152, 1, 6, ncols, beta);

                    simdgeo::geom_i_y(buffer, 16652, 9548, 12152, 1, 6, ncols, beta);

                    simdgeo::geom_i_z(buffer, 16820, 9548, 12152, 1, 6, ncols, beta);

                    simdgeo::geom_i_x(buffer, 16988, 10178, 12800, 1, 6, ncols, beta);

                    simdgeo::geom_i_y(buffer, 17156, 10178, 12800, 1, 6, ncols, beta);

                    simdgeo::geom_i_z(buffer, 17324, 10178, 12800, 1, 6, ncols, beta);

                    simdgeo::geom_k_x(buffer, 17492, 10808, 13448, 1, 6, ncols, beta);

                    simdgeo::geom_k_y(buffer, 17708, 10808, 13448, 1, 6, ncols, beta);

                    simdgeo::geom_k_z(buffer, 17924, 10808, 13448, 1, 6, ncols, beta);

                    simdgeo::geom_k_x(buffer, 18140, 11480, 13988, 1, 6, ncols, beta);

                    simdgeo::geom_k_y(buffer, 18356, 11480, 13988, 1, 6, ncols, beta);

                    simdgeo::geom_k_z(buffer, 18572, 11480, 13988, 1, 6, ncols, beta);

                    simdgeo::geom_l_x(buffer, 18788, 12152, 14528, 1, 6, ncols, beta);

                    simdgeo::geom_l_y(buffer, 19058, 12152, 14528, 1, 6, ncols, beta);

                    simdgeo::geom_l_z(buffer, 19328, 12152, 14528, 1, 6, ncols, beta);

                    simdgeo::geom_l_x(buffer, 19598, 12800, 14858, 1, 6, ncols, beta);

                    simdgeo::geom_l_y(buffer, 19868, 12800, 14858, 1, 6, ncols, beta);

                    simdgeo::geom_l_z(buffer, 20138, 12800, 14858, 1, 6, ncols, beta);

                    simdfunc::contract_primitives(buffer, 20408, 15188, 90, ncols);

                    simdfunc::contract_primitives(buffer, 20573, 15278, 90, ncols);

                    simdfunc::contract_primitives(buffer, 20738, 15368, 90, ncols);

                    simdfunc::contract_primitives(buffer, 20903, 8468, 90, ncols);

                    simdfunc::contract_primitives(buffer, 21068, 15458, 90, ncols);

                    simdfunc::contract_primitives(buffer, 21233, 15548, 90, ncols);

                    simdfunc::contract_primitives(buffer, 21398, 15638, 90, ncols);

                    simdfunc::contract_primitives(buffer, 21563, 9008, 90, ncols);

                    simdfunc::contract_primitives(buffer, 21728, 15728, 126, ncols);

                    simdfunc::contract_primitives(buffer, 21959, 15854, 126, ncols);

                    simdfunc::contract_primitives(buffer, 22190, 15980, 126, ncols);

                    simdfunc::contract_primitives(buffer, 22421, 9548, 126, ncols);

                    simdfunc::contract_primitives(buffer, 22652, 16106, 126, ncols);

                    simdfunc::contract_primitives(buffer, 22883, 16232, 126, ncols);

                    simdfunc::contract_primitives(buffer, 23114, 16358, 126, ncols);

                    simdfunc::contract_primitives(buffer, 23345, 10178, 126, ncols);

                    simdfunc::contract_primitives(buffer, 23576, 16484, 168, ncols);

                    simdfunc::contract_primitives(buffer, 23884, 16652, 168, ncols);

                    simdfunc::contract_primitives(buffer, 24192, 16820, 168, ncols);

                    simdfunc::contract_primitives(buffer, 24500, 10808, 168, ncols);

                    simdfunc::contract_primitives(buffer, 24808, 16988, 168, ncols);

                    simdfunc::contract_primitives(buffer, 25116, 17156, 168, ncols);

                    simdfunc::contract_primitives(buffer, 25424, 17324, 168, ncols);

                    simdfunc::contract_primitives(buffer, 25732, 11480, 168, ncols);

                    simdfunc::contract_primitives(buffer, 26040, 17492, 216, ncols);

                    simdfunc::contract_primitives(buffer, 26436, 17708, 216, ncols);

                    simdfunc::contract_primitives(buffer, 26832, 17924, 216, ncols);

                    simdfunc::contract_primitives(buffer, 27228, 12152, 216, ncols);

                    simdfunc::contract_primitives(buffer, 27624, 18140, 216, ncols);

                    simdfunc::contract_primitives(buffer, 28020, 18356, 216, ncols);

                    simdfunc::contract_primitives(buffer, 28416, 18572, 216, ncols);

                    simdfunc::contract_primitives(buffer, 28812, 12800, 216, ncols);

                    simdfunc::contract_primitives(buffer, 29208, 18788, 270, ncols);

                    simdfunc::contract_primitives(buffer, 29703, 19058, 270, ncols);

                    simdfunc::contract_primitives(buffer, 30198, 19328, 270, ncols);

                    simdfunc::contract_primitives(buffer, 30693, 19598, 270, ncols);

                    simdfunc::contract_primitives(buffer, 31188, 19868, 270, ncols);

                    simdfunc::contract_primitives(buffer, 31683, 20138, 270, ncols);
                }
            }
        }

        simdtrf::transform_d_inner(buffer, 20498, 20408, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 20663, 20573, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 20828, 20738, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 20993, 20903, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 21158, 21068, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 21323, 21233, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 21488, 21398, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 21653, 21563, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 21854, 21728, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 22085, 21959, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 22316, 22190, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 22547, 22421, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 22778, 22652, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 23009, 22883, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 23240, 23114, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 23471, 23345, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 23744, 23576, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 24052, 23884, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 24360, 24192, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 24668, 24500, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 24976, 24808, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 25284, 25116, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 25592, 25424, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 25900, 25732, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 26256, 26040, 36, 1, nmax);

        simdtrf::transform_d_inner(buffer, 26652, 26436, 36, 1, nmax);

        simdtrf::transform_d_inner(buffer, 27048, 26832, 36, 1, nmax);

        simdtrf::transform_d_inner(buffer, 27444, 27228, 36, 1, nmax);

        simdtrf::transform_d_inner(buffer, 27840, 27624, 36, 1, nmax);

        simdtrf::transform_d_inner(buffer, 28236, 28020, 36, 1, nmax);

        simdtrf::transform_d_inner(buffer, 28632, 28416, 36, 1, nmax);

        simdtrf::transform_d_inner(buffer, 29028, 28812, 36, 1, nmax);

        simdtrf::transform_d_inner(buffer, 29478, 29208, 45, 1, nmax);

        simdtrf::transform_d_inner(buffer, 29973, 29703, 45, 1, nmax);

        simdtrf::transform_d_inner(buffer, 30468, 30198, 45, 1, nmax);

        simdtrf::transform_d_inner(buffer, 30963, 30693, 45, 1, nmax);

        simdtrf::transform_d_inner(buffer, 31458, 31188, 45, 1, nmax);

        simdtrf::transform_d_inner(buffer, 31953, 31683, 45, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 32178, 20498, 20993, 21854, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 32403, 20663, 20993, 22085, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 32628, 20828, 20993, 22316, 5,
                                          nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 32853, 20993, 22547, 5, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 33078, 21158, 21653, 22778, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 33303, 21323, 21653, 23009, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 33528, 21488, 21653, 23240, 5,
                                          nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 33753, 21653, 23471, 5, nmax);

        simdtrf::compute_hrr_geom_010x_ph(buffer, coordinates, 33978, 21854, 22547, 23744, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_ph(buffer, coordinates, 34293, 22085, 22547, 24052, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_ph(buffer, coordinates, 34608, 22316, 22547, 24360, 5,
                                          nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 34923, 22547, 24668, 5, nmax);

        simdtrf::compute_hrr_geom_010x_ph(buffer, coordinates, 35238, 22778, 23471, 24976, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_ph(buffer, coordinates, 35553, 23009, 23471, 25284, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_ph(buffer, coordinates, 35868, 23240, 23471, 25592, 5,
                                          nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 36183, 23471, 25900, 5, nmax);

        simdtrf::compute_hrr_geom_010x_pi(buffer, coordinates, 36498, 23744, 24668, 26256, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pi(buffer, coordinates, 36918, 24052, 24668, 26652, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pi(buffer, coordinates, 37338, 24360, 24668, 27048, 5,
                                          nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 37758, 24668, 27444, 5, nmax);

        simdtrf::compute_hrr_geom_010x_pi(buffer, coordinates, 38178, 24976, 25900, 27840, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pi(buffer, coordinates, 38598, 25284, 25900, 28236, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pi(buffer, coordinates, 39018, 25592, 25900, 28632, 5,
                                          nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 39438, 25900, 29028, 5, nmax);

        simdtrf::compute_hrr_geom_010x_pk(buffer, coordinates, 39858, 26256, 27444, 29478, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pk(buffer, coordinates, 40398, 26652, 27444, 29973, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pk(buffer, coordinates, 40938, 27048, 27444, 30468, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_pk(buffer, coordinates, 41478, 27840, 29028, 30963, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pk(buffer, coordinates, 42018, 28236, 29028, 31458, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pk(buffer, coordinates, 42558, 28632, 29028, 31953, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_dg(buffer, coordinates, 43098, 32178, 32853, 33978, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dg(buffer, coordinates, 43548, 32403, 32853, 34293, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dg(buffer, coordinates, 43998, 32628, 32853, 34608, 5,
                                          nmax);

        simdtrf::compute_hrr_dg(buffer, coordinates, 44448, 32853, 34923, 5, nmax);

        simdtrf::compute_hrr_geom_010x_dg(buffer, coordinates, 44898, 33078, 33753, 35238, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dg(buffer, coordinates, 45348, 33303, 33753, 35553, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dg(buffer, coordinates, 45798, 33528, 33753, 35868, 5,
                                          nmax);

        simdtrf::compute_hrr_dg(buffer, coordinates, 46248, 33753, 36183, 5, nmax);

        simdtrf::compute_hrr_geom_010x_dh(buffer, coordinates, 46698, 33978, 34923, 36498, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dh(buffer, coordinates, 47328, 34293, 34923, 36918, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dh(buffer, coordinates, 47958, 34608, 34923, 37338, 5,
                                          nmax);

        simdtrf::compute_hrr_dh(buffer, coordinates, 48588, 34923, 37758, 5, nmax);

        simdtrf::compute_hrr_geom_010x_dh(buffer, coordinates, 49218, 35238, 36183, 38178, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dh(buffer, coordinates, 49848, 35553, 36183, 38598, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dh(buffer, coordinates, 50478, 35868, 36183, 39018, 5,
                                          nmax);

        simdtrf::compute_hrr_dh(buffer, coordinates, 51108, 36183, 39438, 5, nmax);

        simdtrf::compute_hrr_geom_010x_di(buffer, coordinates, 51738, 36498, 37758, 39858, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_di(buffer, coordinates, 52578, 36918, 37758, 40398, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_di(buffer, coordinates, 53418, 37338, 37758, 40938, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_di(buffer, coordinates, 54258, 38178, 39438, 41478, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_di(buffer, coordinates, 55098, 38598, 39438, 42018, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_di(buffer, coordinates, 55938, 39018, 39438, 42558, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_fg(buffer, coordinates, 56778, 43098, 44448, 46698, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_fg(buffer, coordinates, 57528, 43548, 44448, 47328, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_fg(buffer, coordinates, 58278, 43998, 44448, 47958, 5,
                                          nmax);

        simdtrf::compute_hrr_fg(buffer, coordinates, 59028, 44448, 48588, 5, nmax);

        simdtrf::compute_hrr_geom_010x_fg(buffer, coordinates, 59778, 44898, 46248, 49218, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_fg(buffer, coordinates, 60528, 45348, 46248, 49848, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_fg(buffer, coordinates, 61278, 45798, 46248, 50478, 5,
                                          nmax);

        simdtrf::compute_hrr_fg(buffer, coordinates, 62028, 46248, 51108, 5, nmax);

        simdtrf::compute_hrr_geom_010x_fh(buffer, coordinates, 62778, 46698, 48588, 51738, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_fh(buffer, coordinates, 63828, 47328, 48588, 52578, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_fh(buffer, coordinates, 64878, 47958, 48588, 53418, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_fh(buffer, coordinates, 65928, 49218, 51108, 54258, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_fh(buffer, coordinates, 66978, 49848, 51108, 55098, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_fh(buffer, coordinates, 68028, 50478, 51108, 55938, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_gg(buffer, coordinates, 69078, 56778, 59028, 62778, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_gg(buffer, coordinates, 70203, 57528, 59028, 63828, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_gg(buffer, coordinates, 71328, 58278, 59028, 64878, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_gg(buffer, coordinates, 72453, 59778, 62028, 65928, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_gg(buffer, coordinates, 73578, 60528, 62028, 66978, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_gg(buffer, coordinates, 74703, 61278, 62028, 68028, 5,
                                          nmax);

        simdtrf::transform_g_inner(buffer, 75828, 72453, 15, 5, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 75828, 45, nmax);

        simdtrf::transform_g_inner(buffer, 75828, 73578, 15, 5, nmax);

        simdtrf::transform_g_outer(values + 405 * nvalues + n * npairs, nvalues, buffer, 75828,
                                   45, nmax);

        simdtrf::transform_g_inner(buffer, 75828, 74703, 15, 5, nmax);

        simdtrf::transform_g_outer(values + 810 * nvalues + n * npairs, nvalues, buffer, 75828,
                                   45, nmax);

        simdtrf::transform_g_inner(buffer, 75828, 69078, 15, 5, nmax);

        simdtrf::transform_g_outer(values + 1215 * nvalues + n * npairs, nvalues, buffer, 75828,
                                   45, nmax);

        simdtrf::transform_g_inner(buffer, 75828, 70203, 15, 5, nmax);

        simdtrf::transform_g_outer(values + 1620 * nvalues + n * npairs, nvalues, buffer, 75828,
                                   45, nmax);

        simdtrf::transform_g_inner(buffer, 75828, 71328, 15, 5, nmax);

        simdtrf::transform_g_outer(values + 2025 * nvalues + n * npairs, nvalues, buffer, 75828,
                                   45, nmax);
    }

    for (size_t m = 0; m < 2430; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
