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


#include "SimdThreeCenterElectronRepulsionGeom010RsRecGPF.hpp"

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
#include "SimdGeometryP1.hpp"
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
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferDD.hpp"
#include "SimdTransferDP.hpp"
#include "SimdTransferFP.hpp"
#include "SimdTransferGeom010XDD.hpp"
#include "SimdTransferGeom010XDF.hpp"
#include "SimdTransferGeom010XDP.hpp"
#include "SimdTransferGeom010XFD.hpp"
#include "SimdTransferGeom010XFP.hpp"
#include "SimdTransferGeom010XGP.hpp"
#include "SimdTransferGeom010XPD.hpp"
#include "SimdTransferGeom010XPF.hpp"
#include "SimdTransferGeom010XPG.hpp"
#include "SimdTransferGeom010XPP.hpp"
#include "SimdTransferGeom010YDD.hpp"
#include "SimdTransferGeom010YDF.hpp"
#include "SimdTransferGeom010YDP.hpp"
#include "SimdTransferGeom010YFD.hpp"
#include "SimdTransferGeom010YFP.hpp"
#include "SimdTransferGeom010YGP.hpp"
#include "SimdTransferGeom010YPD.hpp"
#include "SimdTransferGeom010YPF.hpp"
#include "SimdTransferGeom010YPG.hpp"
#include "SimdTransferGeom010YPP.hpp"
#include "SimdTransferGeom010ZDD.hpp"
#include "SimdTransferGeom010ZDF.hpp"
#include "SimdTransferGeom010ZDP.hpp"
#include "SimdTransferGeom010ZFD.hpp"
#include "SimdTransferGeom010ZFP.hpp"
#include "SimdTransferGeom010ZGP.hpp"
#include "SimdTransferGeom010ZPD.hpp"
#include "SimdTransferGeom010ZPF.hpp"
#include "SimdTransferGeom010ZPG.hpp"
#include "SimdTransferGeom010ZPP.hpp"
#include "SimdTransferPD.hpp"
#include "SimdTransferPF.hpp"
#include "SimdTransferPP.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_010_gpf_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_010_gpf_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 36813, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1134 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 36813, 13016, 6619, dimensions);

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
                                                            4, 5, 6, 7, 8, 9}, ncols, fj,
                                                            i * nprim_b + j, fq, omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 16, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9}, ncols, fj,
                                                        i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 26, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 29, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 32, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 35, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 38, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 41, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 44, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 47, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 50, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 53, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 56, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 59, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 62, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 65, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 68, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 71, 0, 3, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 74, 0, 3, 7, 8,
                                                                       26, 29, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 80, 0, 3, 8, 9,
                                                                       29, 32, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 86, 0, 3, 9, 10,
                                                                       32, 35, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 92, 0, 3, 10, 11,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 98, 0, 3, 11, 12,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 104, 0, 3, 12, 13,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 110, 0, 3, 13, 14,
                                                                       44, 47, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 116, 0, 3, 17, 18,
                                                                       50, 53, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 122, 0, 3, 18, 19,
                                                                       53, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 128, 0, 3, 19, 20,
                                                                       56, 59, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 134, 0, 3, 20, 21,
                                                                       59, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 140, 0, 3, 21, 22,
                                                                       62, 65, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 146, 0, 3, 22, 23,
                                                                       65, 68, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 152, 0, 3, 23, 24,
                                                                       68, 71, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 158, 0, 3, 26, 29,
                                                                       74, 80, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 168, 0, 3, 29, 32,
                                                                       80, 86, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 178, 0, 3, 32, 35,
                                                                       86, 92, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 188, 0, 3, 35, 38,
                                                                       92, 98, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 198, 0, 3, 38, 41,
                                                                       98, 104, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 208, 0, 3, 41, 44,
                                                                       104, 110, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 218, 0, 3, 50, 53,
                                                                       116, 122, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 228, 0, 3, 53, 56,
                                                                       122, 128, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 238, 0, 3, 56, 59,
                                                                       128, 134, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 248, 0, 3, 59, 62,
                                                                       134, 140, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 258, 0, 3, 62, 65,
                                                                       140, 146, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 268, 0, 3, 65, 68,
                                                                       146, 152, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 278, 0, 3, 74, 80,
                                                                       158, 168, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 293, 0, 3, 80, 86,
                                                                       168, 178, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 308, 0, 3, 86, 92,
                                                                       178, 188, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 323, 0, 3, 92, 98,
                                                                       188, 198, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 338, 0, 3, 98,
                                                                       104, 198, 208, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 353, 0, 3, 116,
                                                                       122, 218, 228, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 368, 0, 3, 122,
                                                                       128, 228, 238, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 383, 0, 3, 128,
                                                                       134, 238, 248, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 398, 0, 3, 134,
                                                                       140, 248, 258, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 413, 0, 3, 140,
                                                                       146, 258, 268, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 428, 0, 3, 158,
                                                                       168, 278, 293, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 449, 0, 3, 168,
                                                                       178, 293, 308, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 470, 0, 3, 178,
                                                                       188, 308, 323, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 491, 0, 3, 188,
                                                                       198, 323, 338, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 512, 0, 3, 218,
                                                                       228, 353, 368, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 533, 0, 3, 228,
                                                                       238, 368, 383, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 554, 0, 3, 238,
                                                                       248, 383, 398, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 575, 0, 3, 248,
                                                                       258, 398, 413, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 596, 0, 3, 278,
                                                                       293, 428, 449, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 624, 0, 3, 293,
                                                                       308, 449, 470, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 652, 0, 3, 308,
                                                                       323, 470, 491, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 680, 0, 3, 353,
                                                                       368, 512, 533, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 708, 0, 3, 368,
                                                                       383, 533, 554, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 736, 0, 3, 383,
                                                                       398, 554, 575, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 764, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 767, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 770, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 773, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 776, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 779, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 782, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 785, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 788, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 791, 3, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 794, 3, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 797, 3, 19, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 800, 3, 20, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 803, 3, 21, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 806, 3, 22, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 809, 3, 23, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 812, 3, 24, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 815, 3, 25, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 818, 3, 9, 32,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 827, 3, 10, 35,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 836, 3, 11, 38,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 845, 3, 12, 41,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 854, 3, 13, 44,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 863, 3, 14, 47,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 872, 3, 19, 56,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 881, 3, 20, 59,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 890, 3, 21, 62,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 899, 3, 22, 65,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 908, 3, 23, 68,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 917, 3, 24, 71,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 926, 3, 26, 74,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 944, 3, 29, 80,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 962, 3, 32, 86,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 980, 3, 35, 92,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 998, 3, 38, 98,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1016, 3, 41, 104,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1034, 3, 44, 110,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1052, 3, 50, 116,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1070, 3, 53, 122,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1088, 3, 56, 128,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1106, 3, 59, 134,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1124, 3, 62, 140,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1142, 3, 65, 146,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1160, 3, 68, 152,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1178, 3, 74, 158,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1208, 3, 80, 168,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1238, 3, 86, 178,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1268, 3, 92, 188,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1298, 3, 98, 198,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1328, 3, 104, 208,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1358, 3, 116, 218,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1388, 3, 122, 228,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1418, 3, 128, 238,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1448, 3, 134, 248,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1478, 3, 140, 258,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1508, 3, 146, 268,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1538, 3, 158, 278,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1583, 3, 168, 293,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1628, 3, 178, 308,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1673, 3, 188, 323,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1718, 3, 198, 338,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1763, 3, 218, 353,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1808, 3, 228, 368,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1853, 3, 238, 383,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1898, 3, 248, 398,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1943, 3, 258, 413,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1988, 3, 278, 428,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2051, 3, 293, 449,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2114, 3, 308, 470,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2177, 3, 323, 491,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2240, 3, 353, 512,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2303, 3, 368, 533,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2366, 3, 383, 554,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2429, 3, 398, 575,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2492, 3, 428, 596,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2576, 3, 449, 624,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2660, 3, 470, 652,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2744, 3, 512, 680,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2828, 3, 533, 708,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2912, 3, 554, 736,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2996, 3, 7, 8,
                                                                       770, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3002, 3, 8, 9,
                                                                       773, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3008, 3, 9, 10,
                                                                       776, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3014, 3, 10, 11,
                                                                       779, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3020, 3, 11, 12,
                                                                       782, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3026, 3, 12, 13,
                                                                       785, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3032, 3, 13, 14,
                                                                       788, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3038, 3, 17, 18,
                                                                       797, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3044, 3, 18, 19,
                                                                       800, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3050, 3, 19, 20,
                                                                       803, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3056, 3, 20, 21,
                                                                       806, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3062, 3, 21, 22,
                                                                       809, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3068, 3, 22, 23,
                                                                       812, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3074, 3, 23, 24,
                                                                       815, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3080, 0, 3, 2996,
                                                                       770, 3002, 818, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3098, 0, 3, 3002,
                                                                       773, 3008, 827, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3116, 0, 3, 3008,
                                                                       776, 3014, 836, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3134, 0, 3, 3014,
                                                                       779, 3020, 845, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3152, 0, 3, 3020,
                                                                       782, 3026, 854, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3170, 0, 3, 3026,
                                                                       785, 3032, 863, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3188, 0, 3, 3038,
                                                                       797, 3044, 872, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3206, 0, 3, 3044,
                                                                       800, 3050, 881, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3224, 0, 3, 3050,
                                                                       803, 3056, 890, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3242, 0, 3, 3056,
                                                                       806, 3062, 899, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3260, 0, 3, 3062,
                                                                       809, 3068, 908, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3278, 0, 3, 3068,
                                                                       812, 3074, 917, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3296, 0, 3, 3080,
                                                                       818, 3098, 74, 80, 962,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3332, 0, 3, 3098,
                                                                       827, 3116, 80, 86, 980,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3368, 0, 3, 3116,
                                                                       836, 3134, 86, 92, 998,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3404, 0, 3, 3134,
                                                                       845, 3152, 92, 98, 1016,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3440, 0, 3, 3152,
                                                                       854, 3170, 98, 104, 1034,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3476, 0, 3, 3188,
                                                                       872, 3206, 116, 122, 1088,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3512, 0, 3, 3206,
                                                                       881, 3224, 122, 128, 1106,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3548, 0, 3, 3224,
                                                                       890, 3242, 128, 134, 1124,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3584, 0, 3, 3242,
                                                                       899, 3260, 134, 140, 1142,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3620, 0, 3, 3260,
                                                                       908, 3278, 140, 146, 1160,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3656, 0, 3, 3296,
                                                                       962, 3332, 158, 168, 1238,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3716, 0, 3, 3332,
                                                                       980, 3368, 168, 178, 1268,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3776, 0, 3, 3368,
                                                                       998, 3404, 178, 188, 1298,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3836, 0, 3, 3404,
                                                                       1016, 3440, 188, 198,
                                                                       1328, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3896, 0, 3, 3476,
                                                                       1088, 3512, 218, 228,
                                                                       1418, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3956, 0, 3, 3512,
                                                                       1106, 3548, 228, 238,
                                                                       1448, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4016, 0, 3, 3548,
                                                                       1124, 3584, 238, 248,
                                                                       1478, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4076, 0, 3, 3584,
                                                                       1142, 3620, 248, 258,
                                                                       1508, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4136, 0, 3, 3656,
                                                                       1238, 3716, 278, 293,
                                                                       1628, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4226, 0, 3, 3716,
                                                                       1268, 3776, 293, 308,
                                                                       1673, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4316, 0, 3, 3776,
                                                                       1298, 3836, 308, 323,
                                                                       1718, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4406, 0, 3, 3896,
                                                                       1418, 3956, 353, 368,
                                                                       1853, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4496, 0, 3, 3956,
                                                                       1448, 4016, 368, 383,
                                                                       1898, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4586, 0, 3, 4016,
                                                                       1478, 4076, 383, 398,
                                                                       1943, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 4676, 0, 3, 4136,
                                                                       1628, 4226, 428, 449,
                                                                       2114, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 4802, 0, 3, 4226,
                                                                       1673, 4316, 449, 470,
                                                                       2177, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 4928, 0, 3, 4406,
                                                                       1853, 4496, 512, 533,
                                                                       2366, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 5054, 0, 3, 4496,
                                                                       1898, 4586, 533, 554,
                                                                       2429, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 5180, 0, 3, 4676,
                                                                       2114, 4802, 596, 624,
                                                                       2660, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 5348, 0, 3, 4928,
                                                                       2366, 5054, 680, 708,
                                                                       2912, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5516, 3, 764, 767,
                                                                       2996, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5526, 3, 767, 770,
                                                                       3002, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5536, 3, 770, 773,
                                                                       3008, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5546, 3, 773, 776,
                                                                       3014, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5556, 3, 776, 779,
                                                                       3020, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5566, 3, 779, 782,
                                                                       3026, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5576, 3, 782, 785,
                                                                       3032, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5586, 3, 791, 794,
                                                                       3038, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5596, 3, 794, 797,
                                                                       3044, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5606, 3, 797, 800,
                                                                       3050, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5616, 3, 800, 803,
                                                                       3056, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5626, 3, 803, 806,
                                                                       3062, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5636, 3, 806, 809,
                                                                       3068, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5646, 3, 809, 812,
                                                                       3074, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5656, 0, 3, 5516,
                                                                       2996, 5526, 3080, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5686, 0, 3, 5526,
                                                                       3002, 5536, 3098, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5716, 0, 3, 5536,
                                                                       3008, 5546, 3116, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5746, 0, 3, 5546,
                                                                       3014, 5556, 3134, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5776, 0, 3, 5556,
                                                                       3020, 5566, 3152, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5806, 0, 3, 5566,
                                                                       3026, 5576, 3170, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5836, 0, 3, 5586,
                                                                       3038, 5596, 3188, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5866, 0, 3, 5596,
                                                                       3044, 5606, 3206, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5896, 0, 3, 5606,
                                                                       3050, 5616, 3224, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5926, 0, 3, 5616,
                                                                       3056, 5626, 3242, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5956, 0, 3, 5626,
                                                                       3062, 5636, 3260, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5986, 0, 3, 5636,
                                                                       3068, 5646, 3278, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 6016, 0, 3, 5656,
                                                                       3080, 5686, 926, 944,
                                                                       3296, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 6076, 0, 3, 5686,
                                                                       3098, 5716, 944, 962,
                                                                       3332, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 6136, 0, 3, 5716,
                                                                       3116, 5746, 962, 980,
                                                                       3368, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 6196, 0, 3, 5746,
                                                                       3134, 5776, 980, 998,
                                                                       3404, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 6256, 0, 3, 5776,
                                                                       3152, 5806, 998, 1016,
                                                                       3440, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 6316, 0, 3, 5836,
                                                                       3188, 5866, 1052, 1070,
                                                                       3476, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 6376, 0, 3, 5866,
                                                                       3206, 5896, 1070, 1088,
                                                                       3512, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 6436, 0, 3, 5896,
                                                                       3224, 5926, 1088, 1106,
                                                                       3548, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 6496, 0, 3, 5926,
                                                                       3242, 5956, 1106, 1124,
                                                                       3584, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 6556, 0, 3, 5956,
                                                                       3260, 5986, 1124, 1142,
                                                                       3620, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 6616, 0, 3, 6016,
                                                                       3296, 6076, 1178, 1208,
                                                                       3656, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 6716, 0, 3, 6076,
                                                                       3332, 6136, 1208, 1238,
                                                                       3716, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 6816, 0, 3, 6136,
                                                                       3368, 6196, 1238, 1268,
                                                                       3776, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 6916, 0, 3, 6196,
                                                                       3404, 6256, 1268, 1298,
                                                                       3836, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 7016, 0, 3, 6316,
                                                                       3476, 6376, 1358, 1388,
                                                                       3896, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 7116, 0, 3, 6376,
                                                                       3512, 6436, 1388, 1418,
                                                                       3956, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 7216, 0, 3, 6436,
                                                                       3548, 6496, 1418, 1448,
                                                                       4016, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 7316, 0, 3, 6496,
                                                                       3584, 6556, 1448, 1478,
                                                                       4076, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 7416, 0, 3, 6616,
                                                                       3656, 6716, 1538, 1583,
                                                                       4136, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 7566, 0, 3, 6716,
                                                                       3716, 6816, 1583, 1628,
                                                                       4226, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 7716, 0, 3, 6816,
                                                                       3776, 6916, 1628, 1673,
                                                                       4316, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 7866, 0, 3, 7016,
                                                                       3896, 7116, 1763, 1808,
                                                                       4406, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 8016, 0, 3, 7116,
                                                                       3956, 7216, 1808, 1853,
                                                                       4496, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 8166, 0, 3, 7216,
                                                                       4016, 7316, 1853, 1898,
                                                                       4586, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 8316, 0, 3, 7416,
                                                                       4136, 7566, 1988, 2051,
                                                                       4676, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 8526, 0, 3, 7566,
                                                                       4226, 7716, 2051, 2114,
                                                                       4802, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 8736, 0, 3, 7866,
                                                                       4406, 8016, 2240, 2303,
                                                                       4928, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 8946, 0, 3, 8016,
                                                                       4496, 8166, 2303, 2366,
                                                                       5054, ncols, gamma, p,
                                                                       q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 9156, 0, 3, 8316,
                                                                       4676, 8526, 2492, 2576,
                                                                       5180, ncols, gamma, p,
                                                                       q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 9436, 0, 3, 8736,
                                                                       4928, 8946, 2744, 2828,
                                                                       5348, ncols, gamma, p,
                                                                       q);

                    simdgeo::geom_p_x(buffer, 9716, 5516, 6016, 1, 10, ncols, beta);

                    simdgeo::geom_p_y(buffer, 9746, 5516, 6016, 1, 10, ncols, beta);

                    simdgeo::geom_p_z(buffer, 9776, 5516, 6016, 1, 10, ncols, beta);

                    simdgeo::geom_p_x(buffer, 9806, 5586, 6316, 1, 10, ncols, beta);

                    simdgeo::geom_p_y(buffer, 9836, 5586, 6316, 1, 10, ncols, beta);

                    simdgeo::geom_p_z(buffer, 9866, 5586, 6316, 1, 10, ncols, beta);

                    simdgeo::geom_d_x(buffer, 9896, 5656, 6616, 1, 10, ncols, beta);

                    simdgeo::geom_d_y(buffer, 9956, 5656, 6616, 1, 10, ncols, beta);

                    simdgeo::geom_d_z(buffer, 10016, 5656, 6616, 1, 10, ncols, beta);

                    simdgeo::geom_d_x(buffer, 10076, 5836, 7016, 1, 10, ncols, beta);

                    simdgeo::geom_d_y(buffer, 10136, 5836, 7016, 1, 10, ncols, beta);

                    simdgeo::geom_d_z(buffer, 10196, 5836, 7016, 1, 10, ncols, beta);

                    simdgeo::geom_f_x(buffer, 10256, 6016, 7416, 1, 10, ncols, beta);

                    simdgeo::geom_f_y(buffer, 10356, 6016, 7416, 1, 10, ncols, beta);

                    simdgeo::geom_f_z(buffer, 10456, 6016, 7416, 1, 10, ncols, beta);

                    simdgeo::geom_f_x(buffer, 10556, 6316, 7866, 1, 10, ncols, beta);

                    simdgeo::geom_f_y(buffer, 10656, 6316, 7866, 1, 10, ncols, beta);

                    simdgeo::geom_f_z(buffer, 10756, 6316, 7866, 1, 10, ncols, beta);

                    simdgeo::geom_g_x(buffer, 10856, 6616, 8316, 1, 10, ncols, beta);

                    simdgeo::geom_g_y(buffer, 11006, 6616, 8316, 1, 10, ncols, beta);

                    simdgeo::geom_g_z(buffer, 11156, 6616, 8316, 1, 10, ncols, beta);

                    simdgeo::geom_g_x(buffer, 11306, 7016, 8736, 1, 10, ncols, beta);

                    simdgeo::geom_g_y(buffer, 11456, 7016, 8736, 1, 10, ncols, beta);

                    simdgeo::geom_g_z(buffer, 11606, 7016, 8736, 1, 10, ncols, beta);

                    simdgeo::geom_h_x(buffer, 11756, 7416, 9156, 1, 10, ncols, beta);

                    simdgeo::geom_h_y(buffer, 11966, 7416, 9156, 1, 10, ncols, beta);

                    simdgeo::geom_h_z(buffer, 12176, 7416, 9156, 1, 10, ncols, beta);

                    simdgeo::geom_h_x(buffer, 12386, 7866, 9436, 1, 10, ncols, beta);

                    simdgeo::geom_h_y(buffer, 12596, 7866, 9436, 1, 10, ncols, beta);

                    simdgeo::geom_h_z(buffer, 12806, 7866, 9436, 1, 10, ncols, beta);

                    simdfunc::contract_primitives(buffer, 13016, 9716, 30, ncols);

                    simdfunc::contract_primitives(buffer, 13067, 9746, 30, ncols);

                    simdfunc::contract_primitives(buffer, 13118, 9776, 30, ncols);

                    simdfunc::contract_primitives(buffer, 13169, 5656, 30, ncols);

                    simdfunc::contract_primitives(buffer, 13220, 9806, 30, ncols);

                    simdfunc::contract_primitives(buffer, 13271, 9836, 30, ncols);

                    simdfunc::contract_primitives(buffer, 13322, 9866, 30, ncols);

                    simdfunc::contract_primitives(buffer, 13373, 5836, 30, ncols);

                    simdfunc::contract_primitives(buffer, 13424, 9896, 60, ncols);

                    simdfunc::contract_primitives(buffer, 13526, 9956, 60, ncols);

                    simdfunc::contract_primitives(buffer, 13628, 10016, 60, ncols);

                    simdfunc::contract_primitives(buffer, 13730, 6016, 60, ncols);

                    simdfunc::contract_primitives(buffer, 13832, 10076, 60, ncols);

                    simdfunc::contract_primitives(buffer, 13934, 10136, 60, ncols);

                    simdfunc::contract_primitives(buffer, 14036, 10196, 60, ncols);

                    simdfunc::contract_primitives(buffer, 14138, 6316, 60, ncols);

                    simdfunc::contract_primitives(buffer, 14240, 10256, 100, ncols);

                    simdfunc::contract_primitives(buffer, 14410, 10356, 100, ncols);

                    simdfunc::contract_primitives(buffer, 14580, 10456, 100, ncols);

                    simdfunc::contract_primitives(buffer, 14750, 6616, 100, ncols);

                    simdfunc::contract_primitives(buffer, 14920, 10556, 100, ncols);

                    simdfunc::contract_primitives(buffer, 15090, 10656, 100, ncols);

                    simdfunc::contract_primitives(buffer, 15260, 10756, 100, ncols);

                    simdfunc::contract_primitives(buffer, 15430, 7016, 100, ncols);

                    simdfunc::contract_primitives(buffer, 15600, 10856, 150, ncols);

                    simdfunc::contract_primitives(buffer, 15855, 11006, 150, ncols);

                    simdfunc::contract_primitives(buffer, 16110, 11156, 150, ncols);

                    simdfunc::contract_primitives(buffer, 16365, 7416, 150, ncols);

                    simdfunc::contract_primitives(buffer, 16620, 11306, 150, ncols);

                    simdfunc::contract_primitives(buffer, 16875, 11456, 150, ncols);

                    simdfunc::contract_primitives(buffer, 17130, 11606, 150, ncols);

                    simdfunc::contract_primitives(buffer, 17385, 7866, 150, ncols);

                    simdfunc::contract_primitives(buffer, 17640, 11756, 210, ncols);

                    simdfunc::contract_primitives(buffer, 17997, 11966, 210, ncols);

                    simdfunc::contract_primitives(buffer, 18354, 12176, 210, ncols);

                    simdfunc::contract_primitives(buffer, 18711, 12386, 210, ncols);

                    simdfunc::contract_primitives(buffer, 19068, 12596, 210, ncols);

                    simdfunc::contract_primitives(buffer, 19425, 12806, 210, ncols);
                }
            }
        }

        simdtrf::transform_f_inner(buffer, 13046, 13016, 3, 1, nmax);

        simdtrf::transform_f_inner(buffer, 13097, 13067, 3, 1, nmax);

        simdtrf::transform_f_inner(buffer, 13148, 13118, 3, 1, nmax);

        simdtrf::transform_f_inner(buffer, 13199, 13169, 3, 1, nmax);

        simdtrf::transform_f_inner(buffer, 13250, 13220, 3, 1, nmax);

        simdtrf::transform_f_inner(buffer, 13301, 13271, 3, 1, nmax);

        simdtrf::transform_f_inner(buffer, 13352, 13322, 3, 1, nmax);

        simdtrf::transform_f_inner(buffer, 13403, 13373, 3, 1, nmax);

        simdtrf::transform_f_inner(buffer, 13484, 13424, 6, 1, nmax);

        simdtrf::transform_f_inner(buffer, 13586, 13526, 6, 1, nmax);

        simdtrf::transform_f_inner(buffer, 13688, 13628, 6, 1, nmax);

        simdtrf::transform_f_inner(buffer, 13790, 13730, 6, 1, nmax);

        simdtrf::transform_f_inner(buffer, 13892, 13832, 6, 1, nmax);

        simdtrf::transform_f_inner(buffer, 13994, 13934, 6, 1, nmax);

        simdtrf::transform_f_inner(buffer, 14096, 14036, 6, 1, nmax);

        simdtrf::transform_f_inner(buffer, 14198, 14138, 6, 1, nmax);

        simdtrf::transform_f_inner(buffer, 14340, 14240, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 14510, 14410, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 14680, 14580, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 14850, 14750, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 15020, 14920, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 15190, 15090, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 15360, 15260, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 15530, 15430, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 15750, 15600, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 16005, 15855, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 16260, 16110, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 16515, 16365, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 16770, 16620, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 17025, 16875, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 17280, 17130, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 17535, 17385, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 17850, 17640, 21, 1, nmax);

        simdtrf::transform_f_inner(buffer, 18207, 17997, 21, 1, nmax);

        simdtrf::transform_f_inner(buffer, 18564, 18354, 21, 1, nmax);

        simdtrf::transform_f_inner(buffer, 18921, 18711, 21, 1, nmax);

        simdtrf::transform_f_inner(buffer, 19278, 19068, 21, 1, nmax);

        simdtrf::transform_f_inner(buffer, 19635, 19425, 21, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pp(buffer, coordinates, 19782, 13046, 13199, 13484, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pp(buffer, coordinates, 19845, 13097, 13199, 13586, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pp(buffer, coordinates, 19908, 13148, 13199, 13688, 7,
                                          nmax);

        simdtrf::compute_hrr_pp(buffer, coordinates, 19971, 13199, 13790, 7, nmax);

        simdtrf::compute_hrr_geom_010x_pp(buffer, coordinates, 20034, 13250, 13403, 13892, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pp(buffer, coordinates, 20097, 13301, 13403, 13994, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pp(buffer, coordinates, 20160, 13352, 13403, 14096, 7,
                                          nmax);

        simdtrf::compute_hrr_pp(buffer, coordinates, 20223, 13403, 14198, 7, nmax);

        simdtrf::compute_hrr_geom_010x_pd(buffer, coordinates, 20286, 13484, 13790, 14340, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pd(buffer, coordinates, 20412, 13586, 13790, 14510, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pd(buffer, coordinates, 20538, 13688, 13790, 14680, 7,
                                          nmax);

        simdtrf::compute_hrr_pd(buffer, coordinates, 20664, 13790, 14850, 7, nmax);

        simdtrf::compute_hrr_geom_010x_pd(buffer, coordinates, 20790, 13892, 14198, 15020, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pd(buffer, coordinates, 20916, 13994, 14198, 15190, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pd(buffer, coordinates, 21042, 14096, 14198, 15360, 7,
                                          nmax);

        simdtrf::compute_hrr_pd(buffer, coordinates, 21168, 14198, 15530, 7, nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 21294, 14340, 14850, 15750, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 21504, 14510, 14850, 16005, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 21714, 14680, 14850, 16260, 7,
                                          nmax);

        simdtrf::compute_hrr_pf(buffer, coordinates, 21924, 14850, 16515, 7, nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 22134, 15020, 15530, 16770, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 22344, 15190, 15530, 17025, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 22554, 15360, 15530, 17280, 7,
                                          nmax);

        simdtrf::compute_hrr_pf(buffer, coordinates, 22764, 15530, 17535, 7, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 22974, 15750, 16515, 17850, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 23289, 16005, 16515, 18207, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 23604, 16260, 16515, 18564, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 23919, 16770, 17535, 18921, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 24234, 17025, 17535, 19278, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 24549, 17280, 17535, 19635, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_dp_out_of_second(buffer, coordinates, 24864, 19782, 19971,
                                                        20286, 7, nmax);

        simdtrf::compute_hrr_geom_010y_dp_out_of_second(buffer, coordinates, 24990, 19845, 19971,
                                                        20412, 7, nmax);

        simdtrf::compute_hrr_geom_010z_dp_out_of_second(buffer, coordinates, 25116, 19908, 19971,
                                                        20538, 7, nmax);

        simdtrf::compute_hrr_dp_out_of_second(buffer, coordinates, 25242, 19971, 20664, 7,
                                              nmax);

        simdtrf::compute_hrr_geom_010x_dp_out_of_second(buffer, coordinates, 25368, 20034, 20223,
                                                        20790, 7, nmax);

        simdtrf::compute_hrr_geom_010y_dp_out_of_second(buffer, coordinates, 25494, 20097, 20223,
                                                        20916, 7, nmax);

        simdtrf::compute_hrr_geom_010z_dp_out_of_second(buffer, coordinates, 25620, 20160, 20223,
                                                        21042, 7, nmax);

        simdtrf::compute_hrr_dp_out_of_second(buffer, coordinates, 25746, 20223, 21168, 7,
                                              nmax);

        simdtrf::compute_hrr_geom_010x_dd(buffer, coordinates, 25872, 20286, 20664, 21294, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dd(buffer, coordinates, 26124, 20412, 20664, 21504, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dd(buffer, coordinates, 26376, 20538, 20664, 21714, 7,
                                          nmax);

        simdtrf::compute_hrr_dd(buffer, coordinates, 26628, 20664, 21924, 7, nmax);

        simdtrf::compute_hrr_geom_010x_dd(buffer, coordinates, 26880, 20790, 21168, 22134, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dd(buffer, coordinates, 27132, 20916, 21168, 22344, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dd(buffer, coordinates, 27384, 21042, 21168, 22554, 7,
                                          nmax);

        simdtrf::compute_hrr_dd(buffer, coordinates, 27636, 21168, 22764, 7, nmax);

        simdtrf::compute_hrr_geom_010x_df(buffer, coordinates, 27888, 21294, 21924, 22974, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_df(buffer, coordinates, 28308, 21504, 21924, 23289, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_df(buffer, coordinates, 28728, 21714, 21924, 23604, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_df(buffer, coordinates, 29148, 22134, 22764, 23919, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_df(buffer, coordinates, 29568, 22344, 22764, 24234, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_df(buffer, coordinates, 29988, 22554, 22764, 24549, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_fp_out_of_second(buffer, coordinates, 30408, 24864, 25242,
                                                        25872, 7, nmax);

        simdtrf::compute_hrr_geom_010y_fp_out_of_second(buffer, coordinates, 30618, 24990, 25242,
                                                        26124, 7, nmax);

        simdtrf::compute_hrr_geom_010z_fp_out_of_second(buffer, coordinates, 30828, 25116, 25242,
                                                        26376, 7, nmax);

        simdtrf::compute_hrr_fp_out_of_second(buffer, coordinates, 31038, 25242, 26628, 7,
                                              nmax);

        simdtrf::compute_hrr_geom_010x_fp_out_of_second(buffer, coordinates, 31248, 25368, 25746,
                                                        26880, 7, nmax);

        simdtrf::compute_hrr_geom_010y_fp_out_of_second(buffer, coordinates, 31458, 25494, 25746,
                                                        27132, 7, nmax);

        simdtrf::compute_hrr_geom_010z_fp_out_of_second(buffer, coordinates, 31668, 25620, 25746,
                                                        27384, 7, nmax);

        simdtrf::compute_hrr_fp_out_of_second(buffer, coordinates, 31878, 25746, 27636, 7,
                                              nmax);

        simdtrf::compute_hrr_geom_010x_fd_out_of_second(buffer, coordinates, 32088, 25872, 26628,
                                                        27888, 7, nmax);

        simdtrf::compute_hrr_geom_010y_fd_out_of_second(buffer, coordinates, 32508, 26124, 26628,
                                                        28308, 7, nmax);

        simdtrf::compute_hrr_geom_010z_fd_out_of_second(buffer, coordinates, 32928, 26376, 26628,
                                                        28728, 7, nmax);

        simdtrf::compute_hrr_geom_010x_fd_out_of_second(buffer, coordinates, 33348, 26880, 27636,
                                                        29148, 7, nmax);

        simdtrf::compute_hrr_geom_010y_fd_out_of_second(buffer, coordinates, 33768, 27132, 27636,
                                                        29568, 7, nmax);

        simdtrf::compute_hrr_geom_010z_fd_out_of_second(buffer, coordinates, 34188, 27384, 27636,
                                                        29988, 7, nmax);

        simdtrf::compute_hrr_geom_010x_gp_out_of_second(buffer, coordinates, 34608, 30408, 31038,
                                                        32088, 7, nmax);

        simdtrf::compute_hrr_geom_010y_gp_out_of_second(buffer, coordinates, 34923, 30618, 31038,
                                                        32508, 7, nmax);

        simdtrf::compute_hrr_geom_010z_gp_out_of_second(buffer, coordinates, 35238, 30828, 31038,
                                                        32928, 7, nmax);

        simdtrf::compute_hrr_geom_010x_gp_out_of_second(buffer, coordinates, 35553, 31248, 31878,
                                                        33348, 7, nmax);

        simdtrf::compute_hrr_geom_010y_gp_out_of_second(buffer, coordinates, 35868, 31458, 31878,
                                                        33768, 7, nmax);

        simdtrf::compute_hrr_geom_010z_gp_out_of_second(buffer, coordinates, 36183, 31668, 31878,
                                                        34188, 7, nmax);

        simdtrf::transform_p_inner(buffer, 36498, 35553, 15, 7, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 36498, 21, nmax);

        simdtrf::transform_p_inner(buffer, 36498, 35868, 15, 7, nmax);

        simdtrf::transform_g_outer(values + 189 * nvalues + n * npairs, nvalues, buffer, 36498,
                                   21, nmax);

        simdtrf::transform_p_inner(buffer, 36498, 36183, 15, 7, nmax);

        simdtrf::transform_g_outer(values + 378 * nvalues + n * npairs, nvalues, buffer, 36498,
                                   21, nmax);

        simdtrf::transform_p_inner(buffer, 36498, 34608, 15, 7, nmax);

        simdtrf::transform_g_outer(values + 567 * nvalues + n * npairs, nvalues, buffer, 36498,
                                   21, nmax);

        simdtrf::transform_p_inner(buffer, 36498, 34923, 15, 7, nmax);

        simdtrf::transform_g_outer(values + 756 * nvalues + n * npairs, nvalues, buffer, 36498,
                                   21, nmax);

        simdtrf::transform_p_inner(buffer, 36498, 35238, 15, 7, nmax);

        simdtrf::transform_g_outer(values + 945 * nvalues + n * npairs, nvalues, buffer, 36498,
                                   21, nmax);
    }

    for (size_t m = 0; m < 1134; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
