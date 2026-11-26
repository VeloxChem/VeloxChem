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


#include "SimdThreeCenterElectronRepulsionGeom010RsRecPDH.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSDD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferGeom010XPD.hpp"
#include "SimdTransferGeom010YPD.hpp"
#include "SimdTransferGeom010ZPD.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_010_pdh_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_010_pdh_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 17997, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 990 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 17997, 13188, 3346, dimensions);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 428, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 431, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 434, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 437, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 440, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 443, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 446, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 449, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 452, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 455, 3, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 458, 3, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 461, 3, 19, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 464, 3, 20, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 467, 3, 21, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 470, 3, 22, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 473, 3, 23, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 476, 3, 24, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 479, 3, 25, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 482, 3, 9, 32,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 491, 3, 10, 35,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 500, 3, 11, 38,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 509, 3, 12, 41,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 518, 3, 13, 44,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 527, 3, 14, 47,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 536, 3, 19, 56,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 545, 3, 20, 59,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 554, 3, 21, 62,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 563, 3, 22, 65,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 572, 3, 23, 68,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 581, 3, 24, 71,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 590, 3, 26, 74,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 608, 3, 29, 80,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 626, 3, 32, 86,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 644, 3, 35, 92,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 662, 3, 38, 98,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 680, 3, 41, 104,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 698, 3, 44, 110,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 716, 3, 50, 116,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 734, 3, 53, 122,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 752, 3, 56, 128,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 770, 3, 59, 134,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 788, 3, 62, 140,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 806, 3, 65, 146,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 824, 3, 68, 152,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 842, 3, 74, 158,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 872, 3, 80, 168,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 902, 3, 86, 178,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 932, 3, 92, 188,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 962, 3, 98, 198,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 992, 3, 104, 208,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1022, 3, 116, 218,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1052, 3, 122, 228,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1082, 3, 128, 238,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1112, 3, 134, 248,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1142, 3, 140, 258,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1172, 3, 146, 268,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1202, 3, 158, 278,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1247, 3, 168, 293,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1292, 3, 178, 308,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1337, 3, 188, 323,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1382, 3, 198, 338,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1427, 3, 218, 353,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1472, 3, 228, 368,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1517, 3, 238, 383,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1562, 3, 248, 398,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1607, 3, 258, 413,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1652, 3, 7, 8,
                                                                       434, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1658, 3, 8, 9,
                                                                       437, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1664, 3, 9, 10,
                                                                       440, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1670, 3, 10, 11,
                                                                       443, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1676, 3, 11, 12,
                                                                       446, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1682, 3, 12, 13,
                                                                       449, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1688, 3, 13, 14,
                                                                       452, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1694, 3, 17, 18,
                                                                       461, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1700, 3, 18, 19,
                                                                       464, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1706, 3, 19, 20,
                                                                       467, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1712, 3, 20, 21,
                                                                       470, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1718, 3, 21, 22,
                                                                       473, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1724, 3, 22, 23,
                                                                       476, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1730, 3, 23, 24,
                                                                       479, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1736, 0, 3, 1652,
                                                                       434, 1658, 482, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1754, 0, 3, 1658,
                                                                       437, 1664, 491, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1772, 0, 3, 1664,
                                                                       440, 1670, 500, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1790, 0, 3, 1670,
                                                                       443, 1676, 509, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1808, 0, 3, 1676,
                                                                       446, 1682, 518, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1826, 0, 3, 1682,
                                                                       449, 1688, 527, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1844, 0, 3, 1694,
                                                                       461, 1700, 536, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1862, 0, 3, 1700,
                                                                       464, 1706, 545, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1880, 0, 3, 1706,
                                                                       467, 1712, 554, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1898, 0, 3, 1712,
                                                                       470, 1718, 563, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1916, 0, 3, 1718,
                                                                       473, 1724, 572, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1934, 0, 3, 1724,
                                                                       476, 1730, 581, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1952, 0, 3, 1736,
                                                                       482, 1754, 74, 80, 626,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1988, 0, 3, 1754,
                                                                       491, 1772, 80, 86, 644,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2024, 0, 3, 1772,
                                                                       500, 1790, 86, 92, 662,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2060, 0, 3, 1790,
                                                                       509, 1808, 92, 98, 680,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2096, 0, 3, 1808,
                                                                       518, 1826, 98, 104, 698,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2132, 0, 3, 1844,
                                                                       536, 1862, 116, 122, 752,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2168, 0, 3, 1862,
                                                                       545, 1880, 122, 128, 770,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2204, 0, 3, 1880,
                                                                       554, 1898, 128, 134, 788,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2240, 0, 3, 1898,
                                                                       563, 1916, 134, 140, 806,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2276, 0, 3, 1916,
                                                                       572, 1934, 140, 146, 824,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2312, 0, 3, 1952,
                                                                       626, 1988, 158, 168, 902,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2372, 0, 3, 1988,
                                                                       644, 2024, 168, 178, 932,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2432, 0, 3, 2024,
                                                                       662, 2060, 178, 188, 962,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2492, 0, 3, 2060,
                                                                       680, 2096, 188, 198, 992,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2552, 0, 3, 2132,
                                                                       752, 2168, 218, 228, 1082,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2612, 0, 3, 2168,
                                                                       770, 2204, 228, 238, 1112,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2672, 0, 3, 2204,
                                                                       788, 2240, 238, 248, 1142,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2732, 0, 3, 2240,
                                                                       806, 2276, 248, 258, 1172,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 2792, 0, 3, 2312,
                                                                       902, 2372, 278, 293, 1292,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 2882, 0, 3, 2372,
                                                                       932, 2432, 293, 308, 1337,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 2972, 0, 3, 2432,
                                                                       962, 2492, 308, 323, 1382,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3062, 0, 3, 2552,
                                                                       1082, 2612, 353, 368,
                                                                       1517, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3152, 0, 3, 2612,
                                                                       1112, 2672, 368, 383,
                                                                       1562, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3242, 0, 3, 2672,
                                                                       1142, 2732, 383, 398,
                                                                       1607, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3332, 3, 428, 431,
                                                                       1652, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3342, 3, 431, 434,
                                                                       1658, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3352, 3, 434, 437,
                                                                       1664, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3362, 3, 437, 440,
                                                                       1670, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3372, 3, 440, 443,
                                                                       1676, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3382, 3, 443, 446,
                                                                       1682, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3392, 3, 446, 449,
                                                                       1688, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3402, 3, 455, 458,
                                                                       1694, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3412, 3, 458, 461,
                                                                       1700, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3422, 3, 461, 464,
                                                                       1706, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3432, 3, 464, 467,
                                                                       1712, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3442, 3, 467, 470,
                                                                       1718, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3452, 3, 470, 473,
                                                                       1724, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3462, 3, 473, 476,
                                                                       1730, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3472, 0, 3, 3332,
                                                                       1652, 3342, 1736, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3502, 0, 3, 3342,
                                                                       1658, 3352, 1754, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3532, 0, 3, 3352,
                                                                       1664, 3362, 1772, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3562, 0, 3, 3362,
                                                                       1670, 3372, 1790, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3592, 0, 3, 3372,
                                                                       1676, 3382, 1808, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3622, 0, 3, 3382,
                                                                       1682, 3392, 1826, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3652, 0, 3, 3402,
                                                                       1694, 3412, 1844, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3682, 0, 3, 3412,
                                                                       1700, 3422, 1862, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3712, 0, 3, 3422,
                                                                       1706, 3432, 1880, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3742, 0, 3, 3432,
                                                                       1712, 3442, 1898, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3772, 0, 3, 3442,
                                                                       1718, 3452, 1916, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3802, 0, 3, 3452,
                                                                       1724, 3462, 1934, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 3832, 0, 3, 3472,
                                                                       1736, 3502, 590, 608,
                                                                       1952, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 3892, 0, 3, 3502,
                                                                       1754, 3532, 608, 626,
                                                                       1988, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 3952, 0, 3, 3532,
                                                                       1772, 3562, 626, 644,
                                                                       2024, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4012, 0, 3, 3562,
                                                                       1790, 3592, 644, 662,
                                                                       2060, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4072, 0, 3, 3592,
                                                                       1808, 3622, 662, 680,
                                                                       2096, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4132, 0, 3, 3652,
                                                                       1844, 3682, 716, 734,
                                                                       2132, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4192, 0, 3, 3682,
                                                                       1862, 3712, 734, 752,
                                                                       2168, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4252, 0, 3, 3712,
                                                                       1880, 3742, 752, 770,
                                                                       2204, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4312, 0, 3, 3742,
                                                                       1898, 3772, 770, 788,
                                                                       2240, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4372, 0, 3, 3772,
                                                                       1916, 3802, 788, 806,
                                                                       2276, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 4432, 0, 3, 3832,
                                                                       1952, 3892, 842, 872,
                                                                       2312, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 4532, 0, 3, 3892,
                                                                       1988, 3952, 872, 902,
                                                                       2372, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 4632, 0, 3, 3952,
                                                                       2024, 4012, 902, 932,
                                                                       2432, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 4732, 0, 3, 4012,
                                                                       2060, 4072, 932, 962,
                                                                       2492, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 4832, 0, 3, 4132,
                                                                       2132, 4192, 1022, 1052,
                                                                       2552, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 4932, 0, 3, 4192,
                                                                       2168, 4252, 1052, 1082,
                                                                       2612, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 5032, 0, 3, 4252,
                                                                       2204, 4312, 1082, 1112,
                                                                       2672, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 5132, 0, 3, 4312,
                                                                       2240, 4372, 1112, 1142,
                                                                       2732, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 5232, 0, 3, 4432,
                                                                       2312, 4532, 1202, 1247,
                                                                       2792, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 5382, 0, 3, 4532,
                                                                       2372, 4632, 1247, 1292,
                                                                       2882, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 5532, 0, 3, 4632,
                                                                       2432, 4732, 1292, 1337,
                                                                       2972, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 5682, 0, 3, 4832,
                                                                       2552, 4932, 1427, 1472,
                                                                       3062, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 5832, 0, 3, 4932,
                                                                       2612, 5032, 1472, 1517,
                                                                       3152, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 5982, 0, 3, 5032,
                                                                       2672, 5132, 1517, 1562,
                                                                       3242, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 6132, 3, 1652,
                                                                       1658, 3352, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 6147, 3, 1658,
                                                                       1664, 3362, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 6162, 3, 1664,
                                                                       1670, 3372, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 6177, 3, 1670,
                                                                       1676, 3382, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 6192, 3, 1676,
                                                                       1682, 3392, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 6207, 3, 1694,
                                                                       1700, 3422, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 6222, 3, 1700,
                                                                       1706, 3432, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 6237, 3, 1706,
                                                                       1712, 3442, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 6252, 3, 1712,
                                                                       1718, 3452, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 6267, 3, 1718,
                                                                       1724, 3462, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 6282, 0, 3, 6132,
                                                                       3352, 6147, 1736, 1754,
                                                                       3532, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 6327, 0, 3, 6147,
                                                                       3362, 6162, 1754, 1772,
                                                                       3562, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 6372, 0, 3, 6162,
                                                                       3372, 6177, 1772, 1790,
                                                                       3592, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 6417, 0, 3, 6177,
                                                                       3382, 6192, 1790, 1808,
                                                                       3622, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 6462, 0, 3, 6207,
                                                                       3422, 6222, 1844, 1862,
                                                                       3712, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 6507, 0, 3, 6222,
                                                                       3432, 6237, 1862, 1880,
                                                                       3742, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 6552, 0, 3, 6237,
                                                                       3442, 6252, 1880, 1898,
                                                                       3772, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 6597, 0, 3, 6252,
                                                                       3452, 6267, 1898, 1916,
                                                                       3802, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 6642, 0, 3, 6282,
                                                                       3532, 6327, 1952, 1988,
                                                                       3952, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 6732, 0, 3, 6327,
                                                                       3562, 6372, 1988, 2024,
                                                                       4012, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 6822, 0, 3, 6372,
                                                                       3592, 6417, 2024, 2060,
                                                                       4072, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 6912, 0, 3, 6462,
                                                                       3712, 6507, 2132, 2168,
                                                                       4252, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 7002, 0, 3, 6507,
                                                                       3742, 6552, 2168, 2204,
                                                                       4312, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 7092, 0, 3, 6552,
                                                                       3772, 6597, 2204, 2240,
                                                                       4372, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 7182, 0, 3, 6642,
                                                                       3952, 6732, 2312, 2372,
                                                                       4632, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 7332, 0, 3, 6732,
                                                                       4012, 6822, 2372, 2432,
                                                                       4732, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 7482, 0, 3, 6912,
                                                                       4252, 7002, 2552, 2612,
                                                                       5032, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 7632, 0, 3, 7002,
                                                                       4312, 7092, 2612, 2672,
                                                                       5132, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 7782, 0, 3, 7182,
                                                                       4632, 7332, 2792, 2882,
                                                                       5532, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 8007, 0, 3, 7482,
                                                                       5032, 7632, 3062, 3152,
                                                                       5982, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 8232, 3, 3332,
                                                                       3342, 6132, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 8253, 3, 3342,
                                                                       3352, 6147, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 8274, 3, 3352,
                                                                       3362, 6162, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 8295, 3, 3362,
                                                                       3372, 6177, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 8316, 3, 3372,
                                                                       3382, 6192, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 8337, 3, 3402,
                                                                       3412, 6207, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 8358, 3, 3412,
                                                                       3422, 6222, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 8379, 3, 3422,
                                                                       3432, 6237, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 8400, 3, 3432,
                                                                       3442, 6252, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 8421, 3, 3442,
                                                                       3452, 6267, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 8442, 0, 3, 8232,
                                                                       6132, 8253, 3472, 3502,
                                                                       6282, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 8505, 0, 3, 8253,
                                                                       6147, 8274, 3502, 3532,
                                                                       6327, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 8568, 0, 3, 8274,
                                                                       6162, 8295, 3532, 3562,
                                                                       6372, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 8631, 0, 3, 8295,
                                                                       6177, 8316, 3562, 3592,
                                                                       6417, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 8694, 0, 3, 8337,
                                                                       6207, 8358, 3652, 3682,
                                                                       6462, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 8757, 0, 3, 8358,
                                                                       6222, 8379, 3682, 3712,
                                                                       6507, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 8820, 0, 3, 8379,
                                                                       6237, 8400, 3712, 3742,
                                                                       6552, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 8883, 0, 3, 8400,
                                                                       6252, 8421, 3742, 3772,
                                                                       6597, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 8946, 0, 3, 8442,
                                                                       6282, 8505, 3832, 3892,
                                                                       6642, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 9072, 0, 3, 8505,
                                                                       6327, 8568, 3892, 3952,
                                                                       6732, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 9198, 0, 3, 8568,
                                                                       6372, 8631, 3952, 4012,
                                                                       6822, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 9324, 0, 3, 8694,
                                                                       6462, 8757, 4132, 4192,
                                                                       6912, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 9450, 0, 3, 8757,
                                                                       6507, 8820, 4192, 4252,
                                                                       7002, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 9576, 0, 3, 8820,
                                                                       6552, 8883, 4252, 4312,
                                                                       7092, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 9702, 0, 3, 8946,
                                                                       6642, 9072, 4432, 4532,
                                                                       7182, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 9912, 0, 3, 9072,
                                                                       6732, 9198, 4532, 4632,
                                                                       7332, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 10122, 0, 3, 9324,
                                                                       6912, 9450, 4832, 4932,
                                                                       7482, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 10332, 0, 3, 9450,
                                                                       7002, 9576, 4932, 5032,
                                                                       7632, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 10542, 0, 3, 9702,
                                                                       7182, 9912, 5232, 5382,
                                                                       7782, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 10857, 0, 3,
                                                                       10122, 7482, 10332, 5682,
                                                                       5832, 8007, ncols, gamma,
                                                                       p, q);

                    simdgeo::geom_d_x(buffer, 11172, 8442, 9702, 1, 21, ncols, beta);

                    simdgeo::geom_d_y(buffer, 11298, 8442, 9702, 1, 21, ncols, beta);

                    simdgeo::geom_d_z(buffer, 11424, 8442, 9702, 1, 21, ncols, beta);

                    simdgeo::geom_d_x(buffer, 11550, 8694, 10122, 1, 21, ncols, beta);

                    simdgeo::geom_d_y(buffer, 11676, 8694, 10122, 1, 21, ncols, beta);

                    simdgeo::geom_d_z(buffer, 11802, 8694, 10122, 1, 21, ncols, beta);

                    simdgeo::geom_f_x(buffer, 11928, 8946, 10542, 1, 21, ncols, beta);

                    simdgeo::geom_f_y(buffer, 12138, 8946, 10542, 1, 21, ncols, beta);

                    simdgeo::geom_f_z(buffer, 12348, 8946, 10542, 1, 21, ncols, beta);

                    simdgeo::geom_f_x(buffer, 12558, 9324, 10857, 1, 21, ncols, beta);

                    simdgeo::geom_f_y(buffer, 12768, 9324, 10857, 1, 21, ncols, beta);

                    simdgeo::geom_f_z(buffer, 12978, 9324, 10857, 1, 21, ncols, beta);

                    simdfunc::contract_primitives(buffer, 13188, 11172, 126, ncols);

                    simdfunc::contract_primitives(buffer, 13380, 11298, 126, ncols);

                    simdfunc::contract_primitives(buffer, 13572, 11424, 126, ncols);

                    simdfunc::contract_primitives(buffer, 13764, 8946, 126, ncols);

                    simdfunc::contract_primitives(buffer, 13956, 11550, 126, ncols);

                    simdfunc::contract_primitives(buffer, 14148, 11676, 126, ncols);

                    simdfunc::contract_primitives(buffer, 14340, 11802, 126, ncols);

                    simdfunc::contract_primitives(buffer, 14532, 9324, 126, ncols);

                    simdfunc::contract_primitives(buffer, 14724, 11928, 210, ncols);

                    simdfunc::contract_primitives(buffer, 15044, 12138, 210, ncols);

                    simdfunc::contract_primitives(buffer, 15364, 12348, 210, ncols);

                    simdfunc::contract_primitives(buffer, 15684, 12558, 210, ncols);

                    simdfunc::contract_primitives(buffer, 16004, 12768, 210, ncols);

                    simdfunc::contract_primitives(buffer, 16324, 12978, 210, ncols);
                }
            }
        }

        simdtrf::transform_h_inner(buffer, 13314, 13188, 6, 1, nmax);

        simdtrf::transform_h_inner(buffer, 13506, 13380, 6, 1, nmax);

        simdtrf::transform_h_inner(buffer, 13698, 13572, 6, 1, nmax);

        simdtrf::transform_h_inner(buffer, 13890, 13764, 6, 1, nmax);

        simdtrf::transform_h_inner(buffer, 14082, 13956, 6, 1, nmax);

        simdtrf::transform_h_inner(buffer, 14274, 14148, 6, 1, nmax);

        simdtrf::transform_h_inner(buffer, 14466, 14340, 6, 1, nmax);

        simdtrf::transform_h_inner(buffer, 14658, 14532, 6, 1, nmax);

        simdtrf::transform_h_inner(buffer, 14934, 14724, 10, 1, nmax);

        simdtrf::transform_h_inner(buffer, 15254, 15044, 10, 1, nmax);

        simdtrf::transform_h_inner(buffer, 15574, 15364, 10, 1, nmax);

        simdtrf::transform_h_inner(buffer, 15894, 15684, 10, 1, nmax);

        simdtrf::transform_h_inner(buffer, 16214, 16004, 10, 1, nmax);

        simdtrf::transform_h_inner(buffer, 16534, 16324, 10, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pd(buffer, coordinates, 16644, 13314, 13890, 14934, 11,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pd(buffer, coordinates, 16842, 13506, 13890, 15254, 11,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pd(buffer, coordinates, 17040, 13698, 13890, 15574, 11,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_pd(buffer, coordinates, 17238, 14082, 14658, 15894, 11,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pd(buffer, coordinates, 17436, 14274, 14658, 16214, 11,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pd(buffer, coordinates, 17634, 14466, 14658, 16534, 11,
                                          nmax);

        simdtrf::transform_d_inner(buffer, 17832, 17238, 3, 11, nmax);

        simdtrf::transform_p_outer(values + n * npairs, nvalues, buffer, 17832, 55, nmax);

        simdtrf::transform_d_inner(buffer, 17832, 17436, 3, 11, nmax);

        simdtrf::transform_p_outer(values + 165 * nvalues + n * npairs, nvalues, buffer, 17832,
                                   55, nmax);

        simdtrf::transform_d_inner(buffer, 17832, 17634, 3, 11, nmax);

        simdtrf::transform_p_outer(values + 330 * nvalues + n * npairs, nvalues, buffer, 17832,
                                   55, nmax);

        simdtrf::transform_d_inner(buffer, 17832, 16644, 3, 11, nmax);

        simdtrf::transform_p_outer(values + 495 * nvalues + n * npairs, nvalues, buffer, 17832,
                                   55, nmax);

        simdtrf::transform_d_inner(buffer, 17832, 16842, 3, 11, nmax);

        simdtrf::transform_p_outer(values + 660 * nvalues + n * npairs, nvalues, buffer, 17832,
                                   55, nmax);

        simdtrf::transform_d_inner(buffer, 17832, 17040, 3, 11, nmax);

        simdtrf::transform_p_outer(values + 825 * nvalues + n * npairs, nvalues, buffer, 17832,
                                   55, nmax);
    }

    for (size_t m = 0; m < 990; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
