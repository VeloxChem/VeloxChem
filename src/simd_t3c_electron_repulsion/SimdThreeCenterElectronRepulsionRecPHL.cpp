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


#include "SimdThreeCenterElectronRepulsionRecPHL.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

#include "SimdThreeCenterElectronRepulsionVrrRecSDD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSID.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSII.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferPH.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformL.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_phl_three_center_electron_repulsion(double               *values,
                                            const size_t          npairs,
                                            const size_t          natoms,
                                            const CBasisFunction &a_function,
                                            const CBasisFunction &b_function,
                                            const CBasisFunction &c_function,
                                            const CSimdMatrix    &coordinates,
                                            const CSimdMatrix    &c_coordinates,
                                            CSimdMatrix          &buffer,
                                            const double          threshold) -> void
{
    if (npairs > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_phl_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto dimensions = simdfunc::make_column_dimensions(
        a_function, b_function, c_function, npairs, coordinates,
        screenfunc::three_center_electron_repulsion_primitive_bound,
        threshold / static_cast<double>(nprims));

    const auto nmax = simdfunc::prepare_buffer(buffer, 61168, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 561 * natoms * npairs, 0.0);

        return;
    }

    const auto pi = mathconst::pi_value();

    // NOTE: a row of the values spans every atom pair of every atom on the ket
    // side, so a kernel handed the block of one atom steps by this to reach the
    // next component -- which is what lets it be the kernel a two-center form
    // uses, unchanged.

    const auto nvalues = natoms * npairs;

    for (size_t n = 0; n < natoms; n++)
    {
        simdfunc::prepare_buffer(buffer, 61168, 56498, 2562, dimensions);

        for (size_t i = 0; i < nprim_a; i++)
        {
            for (size_t j = 0; j < nprim_b; j++)
            {
                const auto p = a_exps[i] + b_exps[j];

                const auto mu = a_exps[i] * b_exps[j] / p;

                const auto fovl = a_norms[i] * b_norms[j];

                const auto fb = a_exps[i] / p;

                const auto fc = b_exps[j] / p;

                simdfunc::compute_pb(buffer, coordinates, 0, nmax, fb);

                simdfunc::compute_pc(buffer, coordinates, c_coordinates, 3, n, nmax, fc);

                simdfunc::compute_pair_exponent(buffer, coordinates, 6, nmax, mu);

                for (size_t k = 0; k < nprim_c; k++)
                {
                    const auto ncols = dimensions[(i * nprim_b + j) * nprim_c + k];

                    if (ncols == 0) continue;

                    const auto gamma = c_exps[k];

                    const auto q = p + gamma;

                    const auto fq = p * gamma / q;

                    const auto fj = 2.0 * fovl * c_norms[k] * pi * pi * std::sqrt(pi)
                                    / (p * gamma * std::sqrt(q));

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 7, 3, 14,
                                                             ncols, fj, 6, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 23, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 26, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 29, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 32, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 35, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 38, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 41, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 44, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 47, 0, 3, 16, 17,
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

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 65, 0, 3, 8, 9,
                                                                       23, 26, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 71, 0, 3, 9, 10,
                                                                       26, 29, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 77, 0, 3, 10, 11,
                                                                       29, 32, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 83, 0, 3, 11, 12,
                                                                       32, 35, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 89, 0, 3, 12, 13,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 95, 0, 3, 13, 14,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 101, 0, 3, 14, 15,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 107, 0, 3, 15, 16,
                                                                       44, 47, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 113, 0, 3, 16, 17,
                                                                       47, 50, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 119, 0, 3, 17, 18,
                                                                       50, 53, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 125, 0, 3, 18, 19,
                                                                       53, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 131, 0, 3, 19, 20,
                                                                       56, 59, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 137, 0, 3, 20, 21,
                                                                       59, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 143, 0, 3, 23, 26,
                                                                       65, 71, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 153, 0, 3, 26, 29,
                                                                       71, 77, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 163, 0, 3, 29, 32,
                                                                       77, 83, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 173, 0, 3, 32, 35,
                                                                       83, 89, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 183, 0, 3, 35, 38,
                                                                       89, 95, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 193, 0, 3, 38, 41,
                                                                       95, 101, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 203, 0, 3, 41, 44,
                                                                       101, 107, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 213, 0, 3, 44, 47,
                                                                       107, 113, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 223, 0, 3, 47, 50,
                                                                       113, 119, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 233, 0, 3, 50, 53,
                                                                       119, 125, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 243, 0, 3, 53, 56,
                                                                       125, 131, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 253, 0, 3, 56, 59,
                                                                       131, 137, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 263, 0, 3, 65, 71,
                                                                       143, 153, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 278, 0, 3, 71, 77,
                                                                       153, 163, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 293, 0, 3, 77, 83,
                                                                       163, 173, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 308, 0, 3, 83, 89,
                                                                       173, 183, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 323, 0, 3, 89, 95,
                                                                       183, 193, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 338, 0, 3, 95,
                                                                       101, 193, 203, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 353, 0, 3, 101,
                                                                       107, 203, 213, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 368, 0, 3, 107,
                                                                       113, 213, 223, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 383, 0, 3, 113,
                                                                       119, 223, 233, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 398, 0, 3, 119,
                                                                       125, 233, 243, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 413, 0, 3, 125,
                                                                       131, 243, 253, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 428, 0, 3, 143,
                                                                       153, 263, 278, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 449, 0, 3, 153,
                                                                       163, 278, 293, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 470, 0, 3, 163,
                                                                       173, 293, 308, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 491, 0, 3, 173,
                                                                       183, 308, 323, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 512, 0, 3, 183,
                                                                       193, 323, 338, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 533, 0, 3, 193,
                                                                       203, 338, 353, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 554, 0, 3, 203,
                                                                       213, 353, 368, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 575, 0, 3, 213,
                                                                       223, 368, 383, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 596, 0, 3, 223,
                                                                       233, 383, 398, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 617, 0, 3, 233,
                                                                       243, 398, 413, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 638, 0, 3, 263,
                                                                       278, 428, 449, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 666, 0, 3, 278,
                                                                       293, 449, 470, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 694, 0, 3, 293,
                                                                       308, 470, 491, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 722, 0, 3, 308,
                                                                       323, 491, 512, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 750, 0, 3, 323,
                                                                       338, 512, 533, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 778, 0, 3, 338,
                                                                       353, 533, 554, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 806, 0, 3, 353,
                                                                       368, 554, 575, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 834, 0, 3, 368,
                                                                       383, 575, 596, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 862, 0, 3, 383,
                                                                       398, 596, 617, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 890, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 893, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 896, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 899, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 902, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 905, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 908, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 911, 3, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 914, 3, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 917, 3, 19, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 920, 3, 20, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 923, 3, 21, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 926, 3, 22, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 929, 3, 10, 29,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 938, 3, 11, 32,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 947, 3, 12, 35,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 956, 3, 13, 38,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 965, 3, 14, 41,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 974, 3, 15, 44,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 983, 3, 16, 47,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 992, 3, 17, 50,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1001, 3, 18, 53,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1010, 3, 19, 56,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1019, 3, 20, 59,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1028, 3, 21, 62,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1037, 3, 29, 77,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1055, 3, 32, 83,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1073, 3, 35, 89,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1091, 3, 38, 95,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1109, 3, 41, 101,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1127, 3, 44, 107,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1145, 3, 47, 113,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1163, 3, 50, 119,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1181, 3, 53, 125,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1199, 3, 56, 131,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1217, 3, 59, 137,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1235, 3, 77, 163,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1265, 3, 83, 173,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1295, 3, 89, 183,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1325, 3, 95, 193,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1355, 3, 101, 203,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1385, 3, 107, 213,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1415, 3, 113, 223,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1445, 3, 119, 233,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1475, 3, 125, 243,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1505, 3, 131, 253,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1535, 3, 163, 293,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1580, 3, 173, 308,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1625, 3, 183, 323,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1670, 3, 193, 338,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1715, 3, 203, 353,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1760, 3, 213, 368,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1805, 3, 223, 383,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1850, 3, 233, 398,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1895, 3, 243, 413,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1940, 3, 293, 470,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2003, 3, 308, 491,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2066, 3, 323, 512,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2129, 3, 338, 533,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2192, 3, 353, 554,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2255, 3, 368, 575,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2318, 3, 383, 596,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2381, 3, 398, 617,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2444, 3, 470, 694,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2528, 3, 491, 722,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2612, 3, 512, 750,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2696, 3, 533, 778,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2780, 3, 554, 806,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2864, 3, 575, 834,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2948, 3, 596, 862,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3032, 3, 8, 9,
                                                                       890, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3038, 3, 9, 10,
                                                                       893, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3044, 3, 10, 11,
                                                                       896, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3050, 3, 11, 12,
                                                                       899, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3056, 3, 12, 13,
                                                                       902, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3062, 3, 13, 14,
                                                                       905, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3068, 3, 14, 15,
                                                                       908, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3074, 3, 15, 16,
                                                                       911, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3080, 3, 16, 17,
                                                                       914, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3086, 3, 17, 18,
                                                                       917, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3092, 3, 18, 19,
                                                                       920, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3098, 3, 19, 20,
                                                                       923, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3104, 3, 20, 21,
                                                                       926, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3110, 0, 3, 3032,
                                                                       890, 3038, 929, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3128, 0, 3, 3038,
                                                                       893, 3044, 938, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3146, 0, 3, 3044,
                                                                       896, 3050, 947, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3164, 0, 3, 3050,
                                                                       899, 3056, 956, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3182, 0, 3, 3056,
                                                                       902, 3062, 965, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3200, 0, 3, 3062,
                                                                       905, 3068, 974, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3218, 0, 3, 3068,
                                                                       908, 3074, 983, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3236, 0, 3, 3074,
                                                                       911, 3080, 992, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3254, 0, 3, 3080,
                                                                       914, 3086, 1001, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3272, 0, 3, 3086,
                                                                       917, 3092, 1010, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3290, 0, 3, 3092,
                                                                       920, 3098, 1019, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3308, 0, 3, 3098,
                                                                       923, 3104, 1028, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3326, 0, 3, 3110,
                                                                       929, 3128, 65, 71, 1037,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3362, 0, 3, 3128,
                                                                       938, 3146, 71, 77, 1055,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3398, 0, 3, 3146,
                                                                       947, 3164, 77, 83, 1073,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3434, 0, 3, 3164,
                                                                       956, 3182, 83, 89, 1091,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3470, 0, 3, 3182,
                                                                       965, 3200, 89, 95, 1109,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3506, 0, 3, 3200,
                                                                       974, 3218, 95, 101, 1127,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3542, 0, 3, 3218,
                                                                       983, 3236, 101, 107, 1145,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3578, 0, 3, 3236,
                                                                       992, 3254, 107, 113, 1163,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3614, 0, 3, 3254,
                                                                       1001, 3272, 113, 119,
                                                                       1181, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3650, 0, 3, 3272,
                                                                       1010, 3290, 119, 125,
                                                                       1199, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3686, 0, 3, 3290,
                                                                       1019, 3308, 125, 131,
                                                                       1217, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3722, 0, 3, 3326,
                                                                       1037, 3362, 143, 153,
                                                                       1235, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3782, 0, 3, 3362,
                                                                       1055, 3398, 153, 163,
                                                                       1265, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3842, 0, 3, 3398,
                                                                       1073, 3434, 163, 173,
                                                                       1295, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3902, 0, 3, 3434,
                                                                       1091, 3470, 173, 183,
                                                                       1325, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3962, 0, 3, 3470,
                                                                       1109, 3506, 183, 193,
                                                                       1355, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4022, 0, 3, 3506,
                                                                       1127, 3542, 193, 203,
                                                                       1385, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4082, 0, 3, 3542,
                                                                       1145, 3578, 203, 213,
                                                                       1415, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4142, 0, 3, 3578,
                                                                       1163, 3614, 213, 223,
                                                                       1445, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4202, 0, 3, 3614,
                                                                       1181, 3650, 223, 233,
                                                                       1475, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4262, 0, 3, 3650,
                                                                       1199, 3686, 233, 243,
                                                                       1505, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4322, 0, 3, 3722,
                                                                       1235, 3782, 263, 278,
                                                                       1535, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4412, 0, 3, 3782,
                                                                       1265, 3842, 278, 293,
                                                                       1580, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4502, 0, 3, 3842,
                                                                       1295, 3902, 293, 308,
                                                                       1625, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4592, 0, 3, 3902,
                                                                       1325, 3962, 308, 323,
                                                                       1670, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4682, 0, 3, 3962,
                                                                       1355, 4022, 323, 338,
                                                                       1715, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4772, 0, 3, 4022,
                                                                       1385, 4082, 338, 353,
                                                                       1760, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4862, 0, 3, 4082,
                                                                       1415, 4142, 353, 368,
                                                                       1805, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4952, 0, 3, 4142,
                                                                       1445, 4202, 368, 383,
                                                                       1850, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5042, 0, 3, 4202,
                                                                       1475, 4262, 383, 398,
                                                                       1895, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 5132, 0, 3, 4322,
                                                                       1535, 4412, 428, 449,
                                                                       1940, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 5258, 0, 3, 4412,
                                                                       1580, 4502, 449, 470,
                                                                       2003, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 5384, 0, 3, 4502,
                                                                       1625, 4592, 470, 491,
                                                                       2066, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 5510, 0, 3, 4592,
                                                                       1670, 4682, 491, 512,
                                                                       2129, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 5636, 0, 3, 4682,
                                                                       1715, 4772, 512, 533,
                                                                       2192, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 5762, 0, 3, 4772,
                                                                       1760, 4862, 533, 554,
                                                                       2255, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 5888, 0, 3, 4862,
                                                                       1805, 4952, 554, 575,
                                                                       2318, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 6014, 0, 3, 4952,
                                                                       1850, 5042, 575, 596,
                                                                       2381, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 6140, 0, 3, 5132,
                                                                       1940, 5258, 638, 666,
                                                                       2444, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 6308, 0, 3, 5258,
                                                                       2003, 5384, 666, 694,
                                                                       2528, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 6476, 0, 3, 5384,
                                                                       2066, 5510, 694, 722,
                                                                       2612, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 6644, 0, 3, 5510,
                                                                       2129, 5636, 722, 750,
                                                                       2696, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 6812, 0, 3, 5636,
                                                                       2192, 5762, 750, 778,
                                                                       2780, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 6980, 0, 3, 5762,
                                                                       2255, 5888, 778, 806,
                                                                       2864, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 7148, 0, 3, 5888,
                                                                       2318, 6014, 806, 834,
                                                                       2948, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7316, 3, 890, 893,
                                                                       3044, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7326, 3, 893, 896,
                                                                       3050, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7336, 3, 896, 899,
                                                                       3056, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7346, 3, 899, 902,
                                                                       3062, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7356, 3, 902, 905,
                                                                       3068, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7366, 3, 905, 908,
                                                                       3074, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7376, 3, 908, 911,
                                                                       3080, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7386, 3, 911, 914,
                                                                       3086, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7396, 3, 914, 917,
                                                                       3092, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7406, 3, 917, 920,
                                                                       3098, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7416, 3, 920, 923,
                                                                       3104, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 7426, 0, 3, 7316,
                                                                       3044, 7326, 3146, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 7456, 0, 3, 7326,
                                                                       3050, 7336, 3164, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 7486, 0, 3, 7336,
                                                                       3056, 7346, 3182, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 7516, 0, 3, 7346,
                                                                       3062, 7356, 3200, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 7546, 0, 3, 7356,
                                                                       3068, 7366, 3218, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 7576, 0, 3, 7366,
                                                                       3074, 7376, 3236, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 7606, 0, 3, 7376,
                                                                       3080, 7386, 3254, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 7636, 0, 3, 7386,
                                                                       3086, 7396, 3272, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 7666, 0, 3, 7396,
                                                                       3092, 7406, 3290, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 7696, 0, 3, 7406,
                                                                       3098, 7416, 3308, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 7726, 0, 3, 7426,
                                                                       3146, 7456, 1037, 1055,
                                                                       3398, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 7786, 0, 3, 7456,
                                                                       3164, 7486, 1055, 1073,
                                                                       3434, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 7846, 0, 3, 7486,
                                                                       3182, 7516, 1073, 1091,
                                                                       3470, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 7906, 0, 3, 7516,
                                                                       3200, 7546, 1091, 1109,
                                                                       3506, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 7966, 0, 3, 7546,
                                                                       3218, 7576, 1109, 1127,
                                                                       3542, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 8026, 0, 3, 7576,
                                                                       3236, 7606, 1127, 1145,
                                                                       3578, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 8086, 0, 3, 7606,
                                                                       3254, 7636, 1145, 1163,
                                                                       3614, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 8146, 0, 3, 7636,
                                                                       3272, 7666, 1163, 1181,
                                                                       3650, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 8206, 0, 3, 7666,
                                                                       3290, 7696, 1181, 1199,
                                                                       3686, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 8266, 0, 3, 7726,
                                                                       3398, 7786, 1235, 1265,
                                                                       3842, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 8366, 0, 3, 7786,
                                                                       3434, 7846, 1265, 1295,
                                                                       3902, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 8466, 0, 3, 7846,
                                                                       3470, 7906, 1295, 1325,
                                                                       3962, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 8566, 0, 3, 7906,
                                                                       3506, 7966, 1325, 1355,
                                                                       4022, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 8666, 0, 3, 7966,
                                                                       3542, 8026, 1355, 1385,
                                                                       4082, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 8766, 0, 3, 8026,
                                                                       3578, 8086, 1385, 1415,
                                                                       4142, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 8866, 0, 3, 8086,
                                                                       3614, 8146, 1415, 1445,
                                                                       4202, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 8966, 0, 3, 8146,
                                                                       3650, 8206, 1445, 1475,
                                                                       4262, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 9066, 0, 3, 8266,
                                                                       3842, 8366, 1535, 1580,
                                                                       4502, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 9216, 0, 3, 8366,
                                                                       3902, 8466, 1580, 1625,
                                                                       4592, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 9366, 0, 3, 8466,
                                                                       3962, 8566, 1625, 1670,
                                                                       4682, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 9516, 0, 3, 8566,
                                                                       4022, 8666, 1670, 1715,
                                                                       4772, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 9666, 0, 3, 8666,
                                                                       4082, 8766, 1715, 1760,
                                                                       4862, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 9816, 0, 3, 8766,
                                                                       4142, 8866, 1760, 1805,
                                                                       4952, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 9966, 0, 3, 8866,
                                                                       4202, 8966, 1805, 1850,
                                                                       5042, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 10116, 0, 3, 9066,
                                                                       4502, 9216, 1940, 2003,
                                                                       5384, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 10326, 0, 3, 9216,
                                                                       4592, 9366, 2003, 2066,
                                                                       5510, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 10536, 0, 3, 9366,
                                                                       4682, 9516, 2066, 2129,
                                                                       5636, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 10746, 0, 3, 9516,
                                                                       4772, 9666, 2129, 2192,
                                                                       5762, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 10956, 0, 3, 9666,
                                                                       4862, 9816, 2192, 2255,
                                                                       5888, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 11166, 0, 3, 9816,
                                                                       4952, 9966, 2255, 2318,
                                                                       6014, ncols, gamma, p,
                                                                       q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 11376, 0, 3,
                                                                       10116, 5384, 10326, 2444,
                                                                       2528, 6476, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 11656, 0, 3,
                                                                       10326, 5510, 10536, 2528,
                                                                       2612, 6644, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 11936, 0, 3,
                                                                       10536, 5636, 10746, 2612,
                                                                       2696, 6812, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 12216, 0, 3,
                                                                       10746, 5762, 10956, 2696,
                                                                       2780, 6980, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 12496, 0, 3,
                                                                       10956, 5888, 11166, 2780,
                                                                       2864, 7148, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 12776, 3, 3032,
                                                                       3038, 7316, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 12791, 3, 3038,
                                                                       3044, 7326, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 12806, 3, 3044,
                                                                       3050, 7336, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 12821, 3, 3050,
                                                                       3056, 7346, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 12836, 3, 3056,
                                                                       3062, 7356, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 12851, 3, 3062,
                                                                       3068, 7366, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 12866, 3, 3068,
                                                                       3074, 7376, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 12881, 3, 3074,
                                                                       3080, 7386, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 12896, 3, 3080,
                                                                       3086, 7396, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 12911, 3, 3086,
                                                                       3092, 7406, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 12926, 3, 3092,
                                                                       3098, 7416, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 12941, 0, 3,
                                                                       12776, 7316, 12791, 3110,
                                                                       3128, 7426, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 12986, 0, 3,
                                                                       12791, 7326, 12806, 3128,
                                                                       3146, 7456, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 13031, 0, 3,
                                                                       12806, 7336, 12821, 3146,
                                                                       3164, 7486, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 13076, 0, 3,
                                                                       12821, 7346, 12836, 3164,
                                                                       3182, 7516, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 13121, 0, 3,
                                                                       12836, 7356, 12851, 3182,
                                                                       3200, 7546, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 13166, 0, 3,
                                                                       12851, 7366, 12866, 3200,
                                                                       3218, 7576, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 13211, 0, 3,
                                                                       12866, 7376, 12881, 3218,
                                                                       3236, 7606, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 13256, 0, 3,
                                                                       12881, 7386, 12896, 3236,
                                                                       3254, 7636, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 13301, 0, 3,
                                                                       12896, 7396, 12911, 3254,
                                                                       3272, 7666, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 13346, 0, 3,
                                                                       12911, 7406, 12926, 3272,
                                                                       3290, 7696, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 13391, 0, 3,
                                                                       12941, 7426, 12986, 3326,
                                                                       3362, 7726, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 13481, 0, 3,
                                                                       12986, 7456, 13031, 3362,
                                                                       3398, 7786, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 13571, 0, 3,
                                                                       13031, 7486, 13076, 3398,
                                                                       3434, 7846, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 13661, 0, 3,
                                                                       13076, 7516, 13121, 3434,
                                                                       3470, 7906, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 13751, 0, 3,
                                                                       13121, 7546, 13166, 3470,
                                                                       3506, 7966, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 13841, 0, 3,
                                                                       13166, 7576, 13211, 3506,
                                                                       3542, 8026, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 13931, 0, 3,
                                                                       13211, 7606, 13256, 3542,
                                                                       3578, 8086, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 14021, 0, 3,
                                                                       13256, 7636, 13301, 3578,
                                                                       3614, 8146, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 14111, 0, 3,
                                                                       13301, 7666, 13346, 3614,
                                                                       3650, 8206, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 14201, 0, 3,
                                                                       13391, 7726, 13481, 3722,
                                                                       3782, 8266, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 14351, 0, 3,
                                                                       13481, 7786, 13571, 3782,
                                                                       3842, 8366, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 14501, 0, 3,
                                                                       13571, 7846, 13661, 3842,
                                                                       3902, 8466, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 14651, 0, 3,
                                                                       13661, 7906, 13751, 3902,
                                                                       3962, 8566, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 14801, 0, 3,
                                                                       13751, 7966, 13841, 3962,
                                                                       4022, 8666, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 14951, 0, 3,
                                                                       13841, 8026, 13931, 4022,
                                                                       4082, 8766, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 15101, 0, 3,
                                                                       13931, 8086, 14021, 4082,
                                                                       4142, 8866, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 15251, 0, 3,
                                                                       14021, 8146, 14111, 4142,
                                                                       4202, 8966, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 15401, 0, 3,
                                                                       14201, 8266, 14351, 4322,
                                                                       4412, 9066, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 15626, 0, 3,
                                                                       14351, 8366, 14501, 4412,
                                                                       4502, 9216, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 15851, 0, 3,
                                                                       14501, 8466, 14651, 4502,
                                                                       4592, 9366, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 16076, 0, 3,
                                                                       14651, 8566, 14801, 4592,
                                                                       4682, 9516, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 16301, 0, 3,
                                                                       14801, 8666, 14951, 4682,
                                                                       4772, 9666, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 16526, 0, 3,
                                                                       14951, 8766, 15101, 4772,
                                                                       4862, 9816, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 16751, 0, 3,
                                                                       15101, 8866, 15251, 4862,
                                                                       4952, 9966, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 16976, 0, 3,
                                                                       15401, 9066, 15626, 5132,
                                                                       5258, 10116, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 17291, 0, 3,
                                                                       15626, 9216, 15851, 5258,
                                                                       5384, 10326, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 17606, 0, 3,
                                                                       15851, 9366, 16076, 5384,
                                                                       5510, 10536, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 17921, 0, 3,
                                                                       16076, 9516, 16301, 5510,
                                                                       5636, 10746, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 18236, 0, 3,
                                                                       16301, 9666, 16526, 5636,
                                                                       5762, 10956, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 18551, 0, 3,
                                                                       16526, 9816, 16751, 5762,
                                                                       5888, 11166, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 18866, 0, 3,
                                                                       16976, 10116, 17291, 6140,
                                                                       6308, 11376, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 19286, 0, 3,
                                                                       17291, 10326, 17606, 6308,
                                                                       6476, 11656, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 19706, 0, 3,
                                                                       17606, 10536, 17921, 6476,
                                                                       6644, 11936, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 20126, 0, 3,
                                                                       17921, 10746, 18236, 6644,
                                                                       6812, 12216, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 20546, 0, 3,
                                                                       18236, 10956, 18551, 6812,
                                                                       6980, 12496, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 20966, 3, 7316,
                                                                       7326, 12806, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 20987, 3, 7326,
                                                                       7336, 12821, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 21008, 3, 7336,
                                                                       7346, 12836, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 21029, 3, 7346,
                                                                       7356, 12851, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 21050, 3, 7356,
                                                                       7366, 12866, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 21071, 3, 7366,
                                                                       7376, 12881, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 21092, 3, 7376,
                                                                       7386, 12896, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 21113, 3, 7386,
                                                                       7396, 12911, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 21134, 3, 7396,
                                                                       7406, 12926, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 21155, 0, 3,
                                                                       20966, 12806, 20987, 7426,
                                                                       7456, 13031, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 21218, 0, 3,
                                                                       20987, 12821, 21008, 7456,
                                                                       7486, 13076, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 21281, 0, 3,
                                                                       21008, 12836, 21029, 7486,
                                                                       7516, 13121, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 21344, 0, 3,
                                                                       21029, 12851, 21050, 7516,
                                                                       7546, 13166, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 21407, 0, 3,
                                                                       21050, 12866, 21071, 7546,
                                                                       7576, 13211, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 21470, 0, 3,
                                                                       21071, 12881, 21092, 7576,
                                                                       7606, 13256, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 21533, 0, 3,
                                                                       21092, 12896, 21113, 7606,
                                                                       7636, 13301, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 21596, 0, 3,
                                                                       21113, 12911, 21134, 7636,
                                                                       7666, 13346, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 21659, 0, 3,
                                                                       21155, 13031, 21218, 7726,
                                                                       7786, 13571, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 21785, 0, 3,
                                                                       21218, 13076, 21281, 7786,
                                                                       7846, 13661, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 21911, 0, 3,
                                                                       21281, 13121, 21344, 7846,
                                                                       7906, 13751, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 22037, 0, 3,
                                                                       21344, 13166, 21407, 7906,
                                                                       7966, 13841, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 22163, 0, 3,
                                                                       21407, 13211, 21470, 7966,
                                                                       8026, 13931, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 22289, 0, 3,
                                                                       21470, 13256, 21533, 8026,
                                                                       8086, 14021, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 22415, 0, 3,
                                                                       21533, 13301, 21596, 8086,
                                                                       8146, 14111, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 22541, 0, 3,
                                                                       21659, 13571, 21785, 8266,
                                                                       8366, 14501, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 22751, 0, 3,
                                                                       21785, 13661, 21911, 8366,
                                                                       8466, 14651, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 22961, 0, 3,
                                                                       21911, 13751, 22037, 8466,
                                                                       8566, 14801, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 23171, 0, 3,
                                                                       22037, 13841, 22163, 8566,
                                                                       8666, 14951, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 23381, 0, 3,
                                                                       22163, 13931, 22289, 8666,
                                                                       8766, 15101, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 23591, 0, 3,
                                                                       22289, 14021, 22415, 8766,
                                                                       8866, 15251, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 23801, 0, 3,
                                                                       22541, 14501, 22751, 9066,
                                                                       9216, 15851, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 24116, 0, 3,
                                                                       22751, 14651, 22961, 9216,
                                                                       9366, 16076, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 24431, 0, 3,
                                                                       22961, 14801, 23171, 9366,
                                                                       9516, 16301, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 24746, 0, 3,
                                                                       23171, 14951, 23381, 9516,
                                                                       9666, 16526, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 25061, 0, 3,
                                                                       23381, 15101, 23591, 9666,
                                                                       9816, 16751, ncols, gamma,
                                                                       p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 25376, 0, 3,
                                                                       23801, 15851, 24116,
                                                                       10116, 10326, 17606,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 25817, 0, 3,
                                                                       24116, 16076, 24431,
                                                                       10326, 10536, 17921,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 26258, 0, 3,
                                                                       24431, 16301, 24746,
                                                                       10536, 10746, 18236,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 26699, 0, 3,
                                                                       24746, 16526, 25061,
                                                                       10746, 10956, 18551,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 27140, 0, 3,
                                                                       25376, 17606, 25817,
                                                                       11376, 11656, 19706,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 27728, 0, 3,
                                                                       25817, 17921, 26258,
                                                                       11656, 11936, 20126,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 28316, 0, 3,
                                                                       26258, 18236, 26699,
                                                                       11936, 12216, 20546,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 28904, 3, 12776,
                                                                       12791, 20966, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 28932, 3, 12791,
                                                                       12806, 20987, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 28960, 3, 12806,
                                                                       12821, 21008, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 28988, 3, 12821,
                                                                       12836, 21029, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 29016, 3, 12836,
                                                                       12851, 21050, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 29044, 3, 12851,
                                                                       12866, 21071, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 29072, 3, 12866,
                                                                       12881, 21092, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 29100, 3, 12881,
                                                                       12896, 21113, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 29128, 3, 12896,
                                                                       12911, 21134, ncols,
                                                                       gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 29156, 0, 3,
                                                                       28904, 20966, 28932,
                                                                       12941, 12986, 21155,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 29240, 0, 3,
                                                                       28932, 20987, 28960,
                                                                       12986, 13031, 21218,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 29324, 0, 3,
                                                                       28960, 21008, 28988,
                                                                       13031, 13076, 21281,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 29408, 0, 3,
                                                                       28988, 21029, 29016,
                                                                       13076, 13121, 21344,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 29492, 0, 3,
                                                                       29016, 21050, 29044,
                                                                       13121, 13166, 21407,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 29576, 0, 3,
                                                                       29044, 21071, 29072,
                                                                       13166, 13211, 21470,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 29660, 0, 3,
                                                                       29072, 21092, 29100,
                                                                       13211, 13256, 21533,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 29744, 0, 3,
                                                                       29100, 21113, 29128,
                                                                       13256, 13301, 21596,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 29828, 0, 3,
                                                                       29156, 21155, 29240,
                                                                       13391, 13481, 21659,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 29996, 0, 3,
                                                                       29240, 21218, 29324,
                                                                       13481, 13571, 21785,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 30164, 0, 3,
                                                                       29324, 21281, 29408,
                                                                       13571, 13661, 21911,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 30332, 0, 3,
                                                                       29408, 21344, 29492,
                                                                       13661, 13751, 22037,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 30500, 0, 3,
                                                                       29492, 21407, 29576,
                                                                       13751, 13841, 22163,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 30668, 0, 3,
                                                                       29576, 21470, 29660,
                                                                       13841, 13931, 22289,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 30836, 0, 3,
                                                                       29660, 21533, 29744,
                                                                       13931, 14021, 22415,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 31004, 0, 3,
                                                                       29828, 21659, 29996,
                                                                       14201, 14351, 22541,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 31284, 0, 3,
                                                                       29996, 21785, 30164,
                                                                       14351, 14501, 22751,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 31564, 0, 3,
                                                                       30164, 21911, 30332,
                                                                       14501, 14651, 22961,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 31844, 0, 3,
                                                                       30332, 22037, 30500,
                                                                       14651, 14801, 23171,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 32124, 0, 3,
                                                                       30500, 22163, 30668,
                                                                       14801, 14951, 23381,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 32404, 0, 3,
                                                                       30668, 22289, 30836,
                                                                       14951, 15101, 23591,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 32684, 0, 3,
                                                                       31004, 22541, 31284,
                                                                       15401, 15626, 23801,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 33104, 0, 3,
                                                                       31284, 22751, 31564,
                                                                       15626, 15851, 24116,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 33524, 0, 3,
                                                                       31564, 22961, 31844,
                                                                       15851, 16076, 24431,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 33944, 0, 3,
                                                                       31844, 23171, 32124,
                                                                       16076, 16301, 24746,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 34364, 0, 3,
                                                                       32124, 23381, 32404,
                                                                       16301, 16526, 25061,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 34784, 0, 3,
                                                                       32684, 23801, 33104,
                                                                       16976, 17291, 25376,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 35372, 0, 3,
                                                                       33104, 24116, 33524,
                                                                       17291, 17606, 25817,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 35960, 0, 3,
                                                                       33524, 24431, 33944,
                                                                       17606, 17921, 26258,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 36548, 0, 3,
                                                                       33944, 24746, 34364,
                                                                       17921, 18236, 26699,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 37136, 0, 3,
                                                                       34784, 25376, 35372,
                                                                       18866, 19286, 27140,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 37920, 0, 3,
                                                                       35372, 25817, 35960,
                                                                       19286, 19706, 27728,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 38704, 0, 3,
                                                                       35960, 26258, 36548,
                                                                       19706, 20126, 28316,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 39488, 3, 20966,
                                                                       20987, 28960, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 39524, 3, 20987,
                                                                       21008, 28988, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 39560, 3, 21008,
                                                                       21029, 29016, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 39596, 3, 21029,
                                                                       21050, 29044, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 39632, 3, 21050,
                                                                       21071, 29072, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 39668, 3, 21071,
                                                                       21092, 29100, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 39704, 3, 21092,
                                                                       21113, 29128, ncols,
                                                                       gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 39740, 0, 3,
                                                                       39488, 28960, 39524,
                                                                       21155, 21218, 29324,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 39848, 0, 3,
                                                                       39524, 28988, 39560,
                                                                       21218, 21281, 29408,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 39956, 0, 3,
                                                                       39560, 29016, 39596,
                                                                       21281, 21344, 29492,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 40064, 0, 3,
                                                                       39596, 29044, 39632,
                                                                       21344, 21407, 29576,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 40172, 0, 3,
                                                                       39632, 29072, 39668,
                                                                       21407, 21470, 29660,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 40280, 0, 3,
                                                                       39668, 29100, 39704,
                                                                       21470, 21533, 29744,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 40388, 0, 3,
                                                                       39740, 29324, 39848,
                                                                       21659, 21785, 30164,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 40604, 0, 3,
                                                                       39848, 29408, 39956,
                                                                       21785, 21911, 30332,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 40820, 0, 3,
                                                                       39956, 29492, 40064,
                                                                       21911, 22037, 30500,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 41036, 0, 3,
                                                                       40064, 29576, 40172,
                                                                       22037, 22163, 30668,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 41252, 0, 3,
                                                                       40172, 29660, 40280,
                                                                       22163, 22289, 30836,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 41468, 0, 3,
                                                                       40388, 30164, 40604,
                                                                       22541, 22751, 31564,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 41828, 0, 3,
                                                                       40604, 30332, 40820,
                                                                       22751, 22961, 31844,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 42188, 0, 3,
                                                                       40820, 30500, 41036,
                                                                       22961, 23171, 32124,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 42548, 0, 3,
                                                                       41036, 30668, 41252,
                                                                       23171, 23381, 32404,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 42908, 0, 3,
                                                                       41468, 31564, 41828,
                                                                       23801, 24116, 33524,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 43448, 0, 3,
                                                                       41828, 31844, 42188,
                                                                       24116, 24431, 33944,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 43988, 0, 3,
                                                                       42188, 32124, 42548,
                                                                       24431, 24746, 34364,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 44528, 0, 3,
                                                                       42908, 33524, 43448,
                                                                       25376, 25817, 35960,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 45284, 0, 3,
                                                                       43448, 33944, 43988,
                                                                       25817, 26258, 36548,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 46040, 0, 3,
                                                                       44528, 35960, 45284,
                                                                       27140, 27728, 38704,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 47048, 3, 28904,
                                                                       28932, 39488, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 47093, 3, 28932,
                                                                       28960, 39524, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 47138, 3, 28960,
                                                                       28988, 39560, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 47183, 3, 28988,
                                                                       29016, 39596, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 47228, 3, 29016,
                                                                       29044, 39632, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 47273, 3, 29044,
                                                                       29072, 39668, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 47318, 3, 29072,
                                                                       29100, 39704, ncols,
                                                                       gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 47363, 0, 3,
                                                                       47048, 39488, 47093,
                                                                       29156, 29240, 39740,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 47498, 0, 3,
                                                                       47093, 39524, 47138,
                                                                       29240, 29324, 39848,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 47633, 0, 3,
                                                                       47138, 39560, 47183,
                                                                       29324, 29408, 39956,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 47768, 0, 3,
                                                                       47183, 39596, 47228,
                                                                       29408, 29492, 40064,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 47903, 0, 3,
                                                                       47228, 39632, 47273,
                                                                       29492, 29576, 40172,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 48038, 0, 3,
                                                                       47273, 39668, 47318,
                                                                       29576, 29660, 40280,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 48173, 0, 3,
                                                                       47363, 39740, 47498,
                                                                       29828, 29996, 40388,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 48443, 0, 3,
                                                                       47498, 39848, 47633,
                                                                       29996, 30164, 40604,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 48713, 0, 3,
                                                                       47633, 39956, 47768,
                                                                       30164, 30332, 40820,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 48983, 0, 3,
                                                                       47768, 40064, 47903,
                                                                       30332, 30500, 41036,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 49253, 0, 3,
                                                                       47903, 40172, 48038,
                                                                       30500, 30668, 41252,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 49523, 0, 3,
                                                                       48173, 40388, 48443,
                                                                       31004, 31284, 41468,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 49973, 0, 3,
                                                                       48443, 40604, 48713,
                                                                       31284, 31564, 41828,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 50423, 0, 3,
                                                                       48713, 40820, 48983,
                                                                       31564, 31844, 42188,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 50873, 0, 3,
                                                                       48983, 41036, 49253,
                                                                       31844, 32124, 42548,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 51323, 0, 3,
                                                                       49523, 41468, 49973,
                                                                       32684, 33104, 42908,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 51998, 0, 3,
                                                                       49973, 41828, 50423,
                                                                       33104, 33524, 43448,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 52673, 0, 3,
                                                                       50423, 42188, 50873,
                                                                       33524, 33944, 43988,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 53348, 0, 3,
                                                                       51323, 42908, 51998,
                                                                       34784, 35372, 44528,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 54293, 0, 3,
                                                                       51998, 43448, 52673,
                                                                       35372, 35960, 45284,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 55238, 0, 3,
                                                                       53348, 44528, 54293,
                                                                       37136, 37920, 46040,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 56498, 53348, 945, ncols);

                    simdfunc::contract_primitives(buffer, 57800, 55238, 1260, ncols);
                }
            }
        }

        simdtrf::transform_l_inner(buffer, 57443, 56498, 21, 1, nmax);

        simdtrf::transform_l_inner(buffer, 59060, 57800, 28, 1, nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 59536, 57443, 59060, 17, nmax);

        simdtrf::transform_h_inner(buffer, 60607, 59536, 3, 17, nmax);

        simdtrf::transform_p_outer(values + n * npairs, nvalues, buffer, 60607, 187, nmax);
    }

    for (size_t m = 0; m < 561; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
