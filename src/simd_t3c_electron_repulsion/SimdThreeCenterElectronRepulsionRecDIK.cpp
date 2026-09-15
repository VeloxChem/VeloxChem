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


#include "SimdThreeCenterElectronRepulsionRecDIK.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSID.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSII.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferDI.hpp"
#include "SimdTransferPI.hpp"
#include "SimdTransferPK.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformI.hpp"
#include "SimdTransformK.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_dik_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_dik_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 97919, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 975 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 97919, 85790, 4884, dimensions);

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

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 7, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15},
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

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 890, 0, 3, 428,
                                                                       449, 638, 666, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 926, 0, 3, 449,
                                                                       470, 666, 694, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 962, 0, 3, 470,
                                                                       491, 694, 722, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 998, 0, 3, 491,
                                                                       512, 722, 750, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1034, 0, 3, 512,
                                                                       533, 750, 778, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1070, 0, 3, 533,
                                                                       554, 778, 806, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1106, 0, 3, 554,
                                                                       575, 806, 834, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1142, 0, 3, 575,
                                                                       596, 834, 862, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1178, 0, 3, 638,
                                                                       666, 890, 926, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1223, 0, 3, 666,
                                                                       694, 926, 962, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1268, 0, 3, 694,
                                                                       722, 962, 998, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1313, 0, 3, 722,
                                                                       750, 998, 1034, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1358, 0, 3, 750,
                                                                       778, 1034, 1070, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1403, 0, 3, 778,
                                                                       806, 1070, 1106, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1448, 0, 3, 806,
                                                                       834, 1106, 1142, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1493, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1496, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1499, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1502, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1505, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1508, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1511, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1514, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1517, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1520, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1523, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1526, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1529, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1532, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1535, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1538, 3, 10, 29,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1547, 3, 11, 32,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1556, 3, 12, 35,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1565, 3, 13, 38,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1574, 3, 14, 41,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1583, 3, 15, 44,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1592, 3, 16, 47,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1601, 3, 17, 50,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1610, 3, 18, 53,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1619, 3, 19, 56,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1628, 3, 20, 59,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1637, 3, 21, 62,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1646, 3, 23, 65,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1664, 3, 26, 71,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1682, 3, 29, 77,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1700, 3, 32, 83,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1718, 3, 35, 89,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1736, 3, 38, 95,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1754, 3, 41, 101,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1772, 3, 44, 107,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1790, 3, 47, 113,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1808, 3, 50, 119,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1826, 3, 53, 125,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1844, 3, 56, 131,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1862, 3, 59, 137,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1880, 3, 65, 143,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1910, 3, 71, 153,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1940, 3, 77, 163,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1970, 3, 83, 173,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2000, 3, 89, 183,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2030, 3, 95, 193,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2060, 3, 101, 203,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2090, 3, 107, 213,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2120, 3, 113, 223,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2150, 3, 119, 233,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2180, 3, 125, 243,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2210, 3, 131, 253,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2240, 3, 143, 263,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2285, 3, 153, 278,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2330, 3, 163, 293,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2375, 3, 173, 308,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2420, 3, 183, 323,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2465, 3, 193, 338,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2510, 3, 203, 353,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2555, 3, 213, 368,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2600, 3, 223, 383,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2645, 3, 233, 398,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2690, 3, 243, 413,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2735, 3, 263, 428,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2798, 3, 278, 449,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2861, 3, 293, 470,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2924, 3, 308, 491,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2987, 3, 323, 512,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3050, 3, 338, 533,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3113, 3, 353, 554,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3176, 3, 368, 575,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3239, 3, 383, 596,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3302, 3, 398, 617,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3365, 3, 428, 638,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3449, 3, 449, 666,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3533, 3, 470, 694,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3617, 3, 491, 722,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3701, 3, 512, 750,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3785, 3, 533, 778,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3869, 3, 554, 806,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3953, 3, 575, 834,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4037, 3, 596, 862,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4121, 3, 638, 890,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4229, 3, 666, 926,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4337, 3, 694, 962,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4445, 3, 722, 998,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4553, 3, 750,
                                                                       1034, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4661, 3, 778,
                                                                       1070, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4769, 3, 806,
                                                                       1106, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4877, 3, 834,
                                                                       1142, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 4985, 3, 890,
                                                                       1178, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5120, 3, 926,
                                                                       1223, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5255, 3, 962,
                                                                       1268, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5390, 3, 998,
                                                                       1313, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5525, 3, 1034,
                                                                       1358, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5660, 3, 1070,
                                                                       1403, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5795, 3, 1106,
                                                                       1448, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5930, 3, 8, 9,
                                                                       1499, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5936, 3, 9, 10,
                                                                       1502, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5942, 3, 10, 11,
                                                                       1505, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5948, 3, 11, 12,
                                                                       1508, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5954, 3, 12, 13,
                                                                       1511, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5960, 3, 13, 14,
                                                                       1514, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5966, 3, 14, 15,
                                                                       1517, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5972, 3, 15, 16,
                                                                       1520, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5978, 3, 16, 17,
                                                                       1523, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5984, 3, 17, 18,
                                                                       1526, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5990, 3, 18, 19,
                                                                       1529, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5996, 3, 19, 20,
                                                                       1532, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6002, 3, 20, 21,
                                                                       1535, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6008, 0, 3, 5930,
                                                                       1499, 5936, 1538, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6026, 0, 3, 5936,
                                                                       1502, 5942, 1547, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6044, 0, 3, 5942,
                                                                       1505, 5948, 1556, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6062, 0, 3, 5948,
                                                                       1508, 5954, 1565, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6080, 0, 3, 5954,
                                                                       1511, 5960, 1574, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6098, 0, 3, 5960,
                                                                       1514, 5966, 1583, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6116, 0, 3, 5966,
                                                                       1517, 5972, 1592, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6134, 0, 3, 5972,
                                                                       1520, 5978, 1601, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6152, 0, 3, 5978,
                                                                       1523, 5984, 1610, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6170, 0, 3, 5984,
                                                                       1526, 5990, 1619, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6188, 0, 3, 5990,
                                                                       1529, 5996, 1628, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6206, 0, 3, 5996,
                                                                       1532, 6002, 1637, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6224, 0, 3, 6008,
                                                                       1538, 6026, 65, 71, 1682,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6260, 0, 3, 6026,
                                                                       1547, 6044, 71, 77, 1700,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6296, 0, 3, 6044,
                                                                       1556, 6062, 77, 83, 1718,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6332, 0, 3, 6062,
                                                                       1565, 6080, 83, 89, 1736,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6368, 0, 3, 6080,
                                                                       1574, 6098, 89, 95, 1754,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6404, 0, 3, 6098,
                                                                       1583, 6116, 95, 101, 1772,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6440, 0, 3, 6116,
                                                                       1592, 6134, 101, 107,
                                                                       1790, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6476, 0, 3, 6134,
                                                                       1601, 6152, 107, 113,
                                                                       1808, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6512, 0, 3, 6152,
                                                                       1610, 6170, 113, 119,
                                                                       1826, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6548, 0, 3, 6170,
                                                                       1619, 6188, 119, 125,
                                                                       1844, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6584, 0, 3, 6188,
                                                                       1628, 6206, 125, 131,
                                                                       1862, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6620, 0, 3, 6224,
                                                                       1682, 6260, 143, 153,
                                                                       1940, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6680, 0, 3, 6260,
                                                                       1700, 6296, 153, 163,
                                                                       1970, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6740, 0, 3, 6296,
                                                                       1718, 6332, 163, 173,
                                                                       2000, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6800, 0, 3, 6332,
                                                                       1736, 6368, 173, 183,
                                                                       2030, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6860, 0, 3, 6368,
                                                                       1754, 6404, 183, 193,
                                                                       2060, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6920, 0, 3, 6404,
                                                                       1772, 6440, 193, 203,
                                                                       2090, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6980, 0, 3, 6440,
                                                                       1790, 6476, 203, 213,
                                                                       2120, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7040, 0, 3, 6476,
                                                                       1808, 6512, 213, 223,
                                                                       2150, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7100, 0, 3, 6512,
                                                                       1826, 6548, 223, 233,
                                                                       2180, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7160, 0, 3, 6548,
                                                                       1844, 6584, 233, 243,
                                                                       2210, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7220, 0, 3, 6620,
                                                                       1940, 6680, 263, 278,
                                                                       2330, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7310, 0, 3, 6680,
                                                                       1970, 6740, 278, 293,
                                                                       2375, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7400, 0, 3, 6740,
                                                                       2000, 6800, 293, 308,
                                                                       2420, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7490, 0, 3, 6800,
                                                                       2030, 6860, 308, 323,
                                                                       2465, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7580, 0, 3, 6860,
                                                                       2060, 6920, 323, 338,
                                                                       2510, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7670, 0, 3, 6920,
                                                                       2090, 6980, 338, 353,
                                                                       2555, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7760, 0, 3, 6980,
                                                                       2120, 7040, 353, 368,
                                                                       2600, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7850, 0, 3, 7040,
                                                                       2150, 7100, 368, 383,
                                                                       2645, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7940, 0, 3, 7100,
                                                                       2180, 7160, 383, 398,
                                                                       2690, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8030, 0, 3, 7220,
                                                                       2330, 7310, 428, 449,
                                                                       2861, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8156, 0, 3, 7310,
                                                                       2375, 7400, 449, 470,
                                                                       2924, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8282, 0, 3, 7400,
                                                                       2420, 7490, 470, 491,
                                                                       2987, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8408, 0, 3, 7490,
                                                                       2465, 7580, 491, 512,
                                                                       3050, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8534, 0, 3, 7580,
                                                                       2510, 7670, 512, 533,
                                                                       3113, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8660, 0, 3, 7670,
                                                                       2555, 7760, 533, 554,
                                                                       3176, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8786, 0, 3, 7760,
                                                                       2600, 7850, 554, 575,
                                                                       3239, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8912, 0, 3, 7850,
                                                                       2645, 7940, 575, 596,
                                                                       3302, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 9038, 0, 3, 8030,
                                                                       2861, 8156, 638, 666,
                                                                       3533, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 9206, 0, 3, 8156,
                                                                       2924, 8282, 666, 694,
                                                                       3617, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 9374, 0, 3, 8282,
                                                                       2987, 8408, 694, 722,
                                                                       3701, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 9542, 0, 3, 8408,
                                                                       3050, 8534, 722, 750,
                                                                       3785, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 9710, 0, 3, 8534,
                                                                       3113, 8660, 750, 778,
                                                                       3869, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 9878, 0, 3, 8660,
                                                                       3176, 8786, 778, 806,
                                                                       3953, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 10046, 0, 3, 8786,
                                                                       3239, 8912, 806, 834,
                                                                       4037, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 10214, 0, 3, 9038,
                                                                       3533, 9206, 890, 926,
                                                                       4337, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 10430, 0, 3, 9206,
                                                                       3617, 9374, 926, 962,
                                                                       4445, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 10646, 0, 3, 9374,
                                                                       3701, 9542, 962, 998,
                                                                       4553, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 10862, 0, 3, 9542,
                                                                       3785, 9710, 998, 1034,
                                                                       4661, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 11078, 0, 3, 9710,
                                                                       3869, 9878, 1034, 1070,
                                                                       4769, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 11294, 0, 3, 9878,
                                                                       3953, 10046, 1070, 1106,
                                                                       4877, ncols, gamma, p,
                                                                       q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 11510, 0, 3,
                                                                       10214, 4337, 10430, 1178,
                                                                       1223, 5255, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 11780, 0, 3,
                                                                       10430, 4445, 10646, 1223,
                                                                       1268, 5390, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 12050, 0, 3,
                                                                       10646, 4553, 10862, 1268,
                                                                       1313, 5525, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 12320, 0, 3,
                                                                       10862, 4661, 11078, 1313,
                                                                       1358, 5660, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 12590, 0, 3,
                                                                       11078, 4769, 11294, 1358,
                                                                       1403, 5795, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12860, 3, 1493,
                                                                       1496, 5930, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12870, 3, 1496,
                                                                       1499, 5936, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12880, 3, 1499,
                                                                       1502, 5942, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12890, 3, 1502,
                                                                       1505, 5948, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12900, 3, 1505,
                                                                       1508, 5954, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12910, 3, 1508,
                                                                       1511, 5960, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12920, 3, 1511,
                                                                       1514, 5966, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12930, 3, 1514,
                                                                       1517, 5972, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12940, 3, 1517,
                                                                       1520, 5978, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12950, 3, 1520,
                                                                       1523, 5984, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12960, 3, 1523,
                                                                       1526, 5990, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12970, 3, 1526,
                                                                       1529, 5996, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12980, 3, 1529,
                                                                       1532, 6002, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12990, 0, 3,
                                                                       12860, 5930, 12870, 6008,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13020, 0, 3,
                                                                       12870, 5936, 12880, 6026,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13050, 0, 3,
                                                                       12880, 5942, 12890, 6044,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13080, 0, 3,
                                                                       12890, 5948, 12900, 6062,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13110, 0, 3,
                                                                       12900, 5954, 12910, 6080,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13140, 0, 3,
                                                                       12910, 5960, 12920, 6098,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13170, 0, 3,
                                                                       12920, 5966, 12930, 6116,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13200, 0, 3,
                                                                       12930, 5972, 12940, 6134,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13230, 0, 3,
                                                                       12940, 5978, 12950, 6152,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13260, 0, 3,
                                                                       12950, 5984, 12960, 6170,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13290, 0, 3,
                                                                       12960, 5990, 12970, 6188,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13320, 0, 3,
                                                                       12970, 5996, 12980, 6206,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13350, 0, 3,
                                                                       12990, 6008, 13020, 1646,
                                                                       1664, 6224, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13410, 0, 3,
                                                                       13020, 6026, 13050, 1664,
                                                                       1682, 6260, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13470, 0, 3,
                                                                       13050, 6044, 13080, 1682,
                                                                       1700, 6296, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13530, 0, 3,
                                                                       13080, 6062, 13110, 1700,
                                                                       1718, 6332, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13590, 0, 3,
                                                                       13110, 6080, 13140, 1718,
                                                                       1736, 6368, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13650, 0, 3,
                                                                       13140, 6098, 13170, 1736,
                                                                       1754, 6404, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13710, 0, 3,
                                                                       13170, 6116, 13200, 1754,
                                                                       1772, 6440, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13770, 0, 3,
                                                                       13200, 6134, 13230, 1772,
                                                                       1790, 6476, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13830, 0, 3,
                                                                       13230, 6152, 13260, 1790,
                                                                       1808, 6512, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13890, 0, 3,
                                                                       13260, 6170, 13290, 1808,
                                                                       1826, 6548, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13950, 0, 3,
                                                                       13290, 6188, 13320, 1826,
                                                                       1844, 6584, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14010, 0, 3,
                                                                       13350, 6224, 13410, 1880,
                                                                       1910, 6620, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14110, 0, 3,
                                                                       13410, 6260, 13470, 1910,
                                                                       1940, 6680, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14210, 0, 3,
                                                                       13470, 6296, 13530, 1940,
                                                                       1970, 6740, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14310, 0, 3,
                                                                       13530, 6332, 13590, 1970,
                                                                       2000, 6800, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14410, 0, 3,
                                                                       13590, 6368, 13650, 2000,
                                                                       2030, 6860, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14510, 0, 3,
                                                                       13650, 6404, 13710, 2030,
                                                                       2060, 6920, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14610, 0, 3,
                                                                       13710, 6440, 13770, 2060,
                                                                       2090, 6980, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14710, 0, 3,
                                                                       13770, 6476, 13830, 2090,
                                                                       2120, 7040, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14810, 0, 3,
                                                                       13830, 6512, 13890, 2120,
                                                                       2150, 7100, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14910, 0, 3,
                                                                       13890, 6548, 13950, 2150,
                                                                       2180, 7160, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 15010, 0, 3,
                                                                       14010, 6620, 14110, 2240,
                                                                       2285, 7220, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 15160, 0, 3,
                                                                       14110, 6680, 14210, 2285,
                                                                       2330, 7310, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 15310, 0, 3,
                                                                       14210, 6740, 14310, 2330,
                                                                       2375, 7400, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 15460, 0, 3,
                                                                       14310, 6800, 14410, 2375,
                                                                       2420, 7490, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 15610, 0, 3,
                                                                       14410, 6860, 14510, 2420,
                                                                       2465, 7580, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 15760, 0, 3,
                                                                       14510, 6920, 14610, 2465,
                                                                       2510, 7670, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 15910, 0, 3,
                                                                       14610, 6980, 14710, 2510,
                                                                       2555, 7760, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 16060, 0, 3,
                                                                       14710, 7040, 14810, 2555,
                                                                       2600, 7850, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 16210, 0, 3,
                                                                       14810, 7100, 14910, 2600,
                                                                       2645, 7940, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 16360, 0, 3,
                                                                       15010, 7220, 15160, 2735,
                                                                       2798, 8030, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 16570, 0, 3,
                                                                       15160, 7310, 15310, 2798,
                                                                       2861, 8156, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 16780, 0, 3,
                                                                       15310, 7400, 15460, 2861,
                                                                       2924, 8282, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 16990, 0, 3,
                                                                       15460, 7490, 15610, 2924,
                                                                       2987, 8408, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 17200, 0, 3,
                                                                       15610, 7580, 15760, 2987,
                                                                       3050, 8534, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 17410, 0, 3,
                                                                       15760, 7670, 15910, 3050,
                                                                       3113, 8660, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 17620, 0, 3,
                                                                       15910, 7760, 16060, 3113,
                                                                       3176, 8786, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 17830, 0, 3,
                                                                       16060, 7850, 16210, 3176,
                                                                       3239, 8912, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 18040, 0, 3,
                                                                       16360, 8030, 16570, 3365,
                                                                       3449, 9038, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 18320, 0, 3,
                                                                       16570, 8156, 16780, 3449,
                                                                       3533, 9206, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 18600, 0, 3,
                                                                       16780, 8282, 16990, 3533,
                                                                       3617, 9374, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 18880, 0, 3,
                                                                       16990, 8408, 17200, 3617,
                                                                       3701, 9542, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 19160, 0, 3,
                                                                       17200, 8534, 17410, 3701,
                                                                       3785, 9710, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 19440, 0, 3,
                                                                       17410, 8660, 17620, 3785,
                                                                       3869, 9878, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 19720, 0, 3,
                                                                       17620, 8786, 17830, 3869,
                                                                       3953, 10046, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 20000, 0, 3,
                                                                       18040, 9038, 18320, 4121,
                                                                       4229, 10214, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 20360, 0, 3,
                                                                       18320, 9206, 18600, 4229,
                                                                       4337, 10430, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 20720, 0, 3,
                                                                       18600, 9374, 18880, 4337,
                                                                       4445, 10646, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 21080, 0, 3,
                                                                       18880, 9542, 19160, 4445,
                                                                       4553, 10862, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 21440, 0, 3,
                                                                       19160, 9710, 19440, 4553,
                                                                       4661, 11078, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 21800, 0, 3,
                                                                       19440, 9878, 19720, 4661,
                                                                       4769, 11294, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 22160, 0, 3,
                                                                       20000, 10214, 20360, 4985,
                                                                       5120, 11510, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 22610, 0, 3,
                                                                       20360, 10430, 20720, 5120,
                                                                       5255, 11780, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 23060, 0, 3,
                                                                       20720, 10646, 21080, 5255,
                                                                       5390, 12050, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 23510, 0, 3,
                                                                       21080, 10862, 21440, 5390,
                                                                       5525, 12320, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 23960, 0, 3,
                                                                       21440, 11078, 21800, 5525,
                                                                       5660, 12590, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24410, 3, 5930,
                                                                       5936, 12880, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24425, 3, 5936,
                                                                       5942, 12890, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24440, 3, 5942,
                                                                       5948, 12900, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24455, 3, 5948,
                                                                       5954, 12910, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24470, 3, 5954,
                                                                       5960, 12920, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24485, 3, 5960,
                                                                       5966, 12930, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24500, 3, 5966,
                                                                       5972, 12940, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24515, 3, 5972,
                                                                       5978, 12950, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24530, 3, 5978,
                                                                       5984, 12960, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24545, 3, 5984,
                                                                       5990, 12970, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24560, 3, 5990,
                                                                       5996, 12980, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 24575, 0, 3,
                                                                       24410, 12880, 24425, 6008,
                                                                       6026, 13050, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 24620, 0, 3,
                                                                       24425, 12890, 24440, 6026,
                                                                       6044, 13080, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 24665, 0, 3,
                                                                       24440, 12900, 24455, 6044,
                                                                       6062, 13110, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 24710, 0, 3,
                                                                       24455, 12910, 24470, 6062,
                                                                       6080, 13140, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 24755, 0, 3,
                                                                       24470, 12920, 24485, 6080,
                                                                       6098, 13170, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 24800, 0, 3,
                                                                       24485, 12930, 24500, 6098,
                                                                       6116, 13200, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 24845, 0, 3,
                                                                       24500, 12940, 24515, 6116,
                                                                       6134, 13230, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 24890, 0, 3,
                                                                       24515, 12950, 24530, 6134,
                                                                       6152, 13260, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 24935, 0, 3,
                                                                       24530, 12960, 24545, 6152,
                                                                       6170, 13290, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 24980, 0, 3,
                                                                       24545, 12970, 24560, 6170,
                                                                       6188, 13320, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 25025, 0, 3,
                                                                       24575, 13050, 24620, 6224,
                                                                       6260, 13470, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 25115, 0, 3,
                                                                       24620, 13080, 24665, 6260,
                                                                       6296, 13530, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 25205, 0, 3,
                                                                       24665, 13110, 24710, 6296,
                                                                       6332, 13590, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 25295, 0, 3,
                                                                       24710, 13140, 24755, 6332,
                                                                       6368, 13650, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 25385, 0, 3,
                                                                       24755, 13170, 24800, 6368,
                                                                       6404, 13710, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 25475, 0, 3,
                                                                       24800, 13200, 24845, 6404,
                                                                       6440, 13770, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 25565, 0, 3,
                                                                       24845, 13230, 24890, 6440,
                                                                       6476, 13830, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 25655, 0, 3,
                                                                       24890, 13260, 24935, 6476,
                                                                       6512, 13890, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 25745, 0, 3,
                                                                       24935, 13290, 24980, 6512,
                                                                       6548, 13950, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 25835, 0, 3,
                                                                       25025, 13470, 25115, 6620,
                                                                       6680, 14210, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 25985, 0, 3,
                                                                       25115, 13530, 25205, 6680,
                                                                       6740, 14310, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 26135, 0, 3,
                                                                       25205, 13590, 25295, 6740,
                                                                       6800, 14410, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 26285, 0, 3,
                                                                       25295, 13650, 25385, 6800,
                                                                       6860, 14510, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 26435, 0, 3,
                                                                       25385, 13710, 25475, 6860,
                                                                       6920, 14610, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 26585, 0, 3,
                                                                       25475, 13770, 25565, 6920,
                                                                       6980, 14710, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 26735, 0, 3,
                                                                       25565, 13830, 25655, 6980,
                                                                       7040, 14810, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 26885, 0, 3,
                                                                       25655, 13890, 25745, 7040,
                                                                       7100, 14910, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 27035, 0, 3,
                                                                       25835, 14210, 25985, 7220,
                                                                       7310, 15310, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 27260, 0, 3,
                                                                       25985, 14310, 26135, 7310,
                                                                       7400, 15460, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 27485, 0, 3,
                                                                       26135, 14410, 26285, 7400,
                                                                       7490, 15610, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 27710, 0, 3,
                                                                       26285, 14510, 26435, 7490,
                                                                       7580, 15760, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 27935, 0, 3,
                                                                       26435, 14610, 26585, 7580,
                                                                       7670, 15910, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 28160, 0, 3,
                                                                       26585, 14710, 26735, 7670,
                                                                       7760, 16060, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 28385, 0, 3,
                                                                       26735, 14810, 26885, 7760,
                                                                       7850, 16210, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 28610, 0, 3,
                                                                       27035, 15310, 27260, 8030,
                                                                       8156, 16780, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 28925, 0, 3,
                                                                       27260, 15460, 27485, 8156,
                                                                       8282, 16990, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 29240, 0, 3,
                                                                       27485, 15610, 27710, 8282,
                                                                       8408, 17200, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 29555, 0, 3,
                                                                       27710, 15760, 27935, 8408,
                                                                       8534, 17410, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 29870, 0, 3,
                                                                       27935, 15910, 28160, 8534,
                                                                       8660, 17620, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 30185, 0, 3,
                                                                       28160, 16060, 28385, 8660,
                                                                       8786, 17830, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 30500, 0, 3,
                                                                       28610, 16780, 28925, 9038,
                                                                       9206, 18600, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 30920, 0, 3,
                                                                       28925, 16990, 29240, 9206,
                                                                       9374, 18880, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 31340, 0, 3,
                                                                       29240, 17200, 29555, 9374,
                                                                       9542, 19160, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 31760, 0, 3,
                                                                       29555, 17410, 29870, 9542,
                                                                       9710, 19440, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 32180, 0, 3,
                                                                       29870, 17620, 30185, 9710,
                                                                       9878, 19720, ncols, gamma,
                                                                       p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 32600, 0, 3,
                                                                       30500, 18600, 30920,
                                                                       10214, 10430, 20720,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 33140, 0, 3,
                                                                       30920, 18880, 31340,
                                                                       10430, 10646, 21080,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 33680, 0, 3,
                                                                       31340, 19160, 31760,
                                                                       10646, 10862, 21440,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 34220, 0, 3,
                                                                       31760, 19440, 32180,
                                                                       10862, 11078, 21800,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 34760, 0, 3,
                                                                       32600, 20720, 33140,
                                                                       11510, 11780, 23060,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 35435, 0, 3,
                                                                       33140, 21080, 33680,
                                                                       11780, 12050, 23510,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 36110, 0, 3,
                                                                       33680, 21440, 34220,
                                                                       12050, 12320, 23960,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 36785, 3, 12860,
                                                                       12870, 24410, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 36806, 3, 12870,
                                                                       12880, 24425, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 36827, 3, 12880,
                                                                       12890, 24440, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 36848, 3, 12890,
                                                                       12900, 24455, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 36869, 3, 12900,
                                                                       12910, 24470, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 36890, 3, 12910,
                                                                       12920, 24485, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 36911, 3, 12920,
                                                                       12930, 24500, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 36932, 3, 12930,
                                                                       12940, 24515, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 36953, 3, 12940,
                                                                       12950, 24530, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 36974, 3, 12950,
                                                                       12960, 24545, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 36995, 3, 12960,
                                                                       12970, 24560, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 37016, 0, 3,
                                                                       36785, 24410, 36806,
                                                                       12990, 13020, 24575,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 37079, 0, 3,
                                                                       36806, 24425, 36827,
                                                                       13020, 13050, 24620,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 37142, 0, 3,
                                                                       36827, 24440, 36848,
                                                                       13050, 13080, 24665,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 37205, 0, 3,
                                                                       36848, 24455, 36869,
                                                                       13080, 13110, 24710,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 37268, 0, 3,
                                                                       36869, 24470, 36890,
                                                                       13110, 13140, 24755,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 37331, 0, 3,
                                                                       36890, 24485, 36911,
                                                                       13140, 13170, 24800,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 37394, 0, 3,
                                                                       36911, 24500, 36932,
                                                                       13170, 13200, 24845,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 37457, 0, 3,
                                                                       36932, 24515, 36953,
                                                                       13200, 13230, 24890,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 37520, 0, 3,
                                                                       36953, 24530, 36974,
                                                                       13230, 13260, 24935,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 37583, 0, 3,
                                                                       36974, 24545, 36995,
                                                                       13260, 13290, 24980,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 37646, 0, 3,
                                                                       37016, 24575, 37079,
                                                                       13350, 13410, 25025,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 37772, 0, 3,
                                                                       37079, 24620, 37142,
                                                                       13410, 13470, 25115,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 37898, 0, 3,
                                                                       37142, 24665, 37205,
                                                                       13470, 13530, 25205,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 38024, 0, 3,
                                                                       37205, 24710, 37268,
                                                                       13530, 13590, 25295,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 38150, 0, 3,
                                                                       37268, 24755, 37331,
                                                                       13590, 13650, 25385,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 38276, 0, 3,
                                                                       37331, 24800, 37394,
                                                                       13650, 13710, 25475,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 38402, 0, 3,
                                                                       37394, 24845, 37457,
                                                                       13710, 13770, 25565,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 38528, 0, 3,
                                                                       37457, 24890, 37520,
                                                                       13770, 13830, 25655,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 38654, 0, 3,
                                                                       37520, 24935, 37583,
                                                                       13830, 13890, 25745,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 38780, 0, 3,
                                                                       37646, 25025, 37772,
                                                                       14010, 14110, 25835,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 38990, 0, 3,
                                                                       37772, 25115, 37898,
                                                                       14110, 14210, 25985,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 39200, 0, 3,
                                                                       37898, 25205, 38024,
                                                                       14210, 14310, 26135,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 39410, 0, 3,
                                                                       38024, 25295, 38150,
                                                                       14310, 14410, 26285,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 39620, 0, 3,
                                                                       38150, 25385, 38276,
                                                                       14410, 14510, 26435,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 39830, 0, 3,
                                                                       38276, 25475, 38402,
                                                                       14510, 14610, 26585,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 40040, 0, 3,
                                                                       38402, 25565, 38528,
                                                                       14610, 14710, 26735,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 40250, 0, 3,
                                                                       38528, 25655, 38654,
                                                                       14710, 14810, 26885,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 40460, 0, 3,
                                                                       38780, 25835, 38990,
                                                                       15010, 15160, 27035,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 40775, 0, 3,
                                                                       38990, 25985, 39200,
                                                                       15160, 15310, 27260,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 41090, 0, 3,
                                                                       39200, 26135, 39410,
                                                                       15310, 15460, 27485,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 41405, 0, 3,
                                                                       39410, 26285, 39620,
                                                                       15460, 15610, 27710,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 41720, 0, 3,
                                                                       39620, 26435, 39830,
                                                                       15610, 15760, 27935,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 42035, 0, 3,
                                                                       39830, 26585, 40040,
                                                                       15760, 15910, 28160,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 42350, 0, 3,
                                                                       40040, 26735, 40250,
                                                                       15910, 16060, 28385,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 42665, 0, 3,
                                                                       40460, 27035, 40775,
                                                                       16360, 16570, 28610,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 43106, 0, 3,
                                                                       40775, 27260, 41090,
                                                                       16570, 16780, 28925,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 43547, 0, 3,
                                                                       41090, 27485, 41405,
                                                                       16780, 16990, 29240,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 43988, 0, 3,
                                                                       41405, 27710, 41720,
                                                                       16990, 17200, 29555,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 44429, 0, 3,
                                                                       41720, 27935, 42035,
                                                                       17200, 17410, 29870,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 44870, 0, 3,
                                                                       42035, 28160, 42350,
                                                                       17410, 17620, 30185,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 45311, 0, 3,
                                                                       42665, 28610, 43106,
                                                                       18040, 18320, 30500,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 45899, 0, 3,
                                                                       43106, 28925, 43547,
                                                                       18320, 18600, 30920,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 46487, 0, 3,
                                                                       43547, 29240, 43988,
                                                                       18600, 18880, 31340,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 47075, 0, 3,
                                                                       43988, 29555, 44429,
                                                                       18880, 19160, 31760,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 47663, 0, 3,
                                                                       44429, 29870, 44870,
                                                                       19160, 19440, 32180,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 48251, 0, 3,
                                                                       45311, 30500, 45899,
                                                                       20000, 20360, 32600,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 49007, 0, 3,
                                                                       45899, 30920, 46487,
                                                                       20360, 20720, 33140,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 49763, 0, 3,
                                                                       46487, 31340, 47075,
                                                                       20720, 21080, 33680,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 50519, 0, 3,
                                                                       47075, 31760, 47663,
                                                                       21080, 21440, 34220,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 51275, 0, 3,
                                                                       48251, 32600, 49007,
                                                                       22160, 22610, 34760,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 52220, 0, 3,
                                                                       49007, 33140, 49763,
                                                                       22610, 23060, 35435,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 53165, 0, 3,
                                                                       49763, 33680, 50519,
                                                                       23060, 23510, 36110,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 54110, 3, 24410,
                                                                       24425, 36827, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 54138, 3, 24425,
                                                                       24440, 36848, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 54166, 3, 24440,
                                                                       24455, 36869, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 54194, 3, 24455,
                                                                       24470, 36890, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 54222, 3, 24470,
                                                                       24485, 36911, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 54250, 3, 24485,
                                                                       24500, 36932, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 54278, 3, 24500,
                                                                       24515, 36953, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 54306, 3, 24515,
                                                                       24530, 36974, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 54334, 3, 24530,
                                                                       24545, 36995, ncols,
                                                                       gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 54362, 0, 3,
                                                                       54110, 36827, 54138,
                                                                       24575, 24620, 37142,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 54446, 0, 3,
                                                                       54138, 36848, 54166,
                                                                       24620, 24665, 37205,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 54530, 0, 3,
                                                                       54166, 36869, 54194,
                                                                       24665, 24710, 37268,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 54614, 0, 3,
                                                                       54194, 36890, 54222,
                                                                       24710, 24755, 37331,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 54698, 0, 3,
                                                                       54222, 36911, 54250,
                                                                       24755, 24800, 37394,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 54782, 0, 3,
                                                                       54250, 36932, 54278,
                                                                       24800, 24845, 37457,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 54866, 0, 3,
                                                                       54278, 36953, 54306,
                                                                       24845, 24890, 37520,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 54950, 0, 3,
                                                                       54306, 36974, 54334,
                                                                       24890, 24935, 37583,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 55034, 0, 3,
                                                                       54362, 37142, 54446,
                                                                       25025, 25115, 37898,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 55202, 0, 3,
                                                                       54446, 37205, 54530,
                                                                       25115, 25205, 38024,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 55370, 0, 3,
                                                                       54530, 37268, 54614,
                                                                       25205, 25295, 38150,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 55538, 0, 3,
                                                                       54614, 37331, 54698,
                                                                       25295, 25385, 38276,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 55706, 0, 3,
                                                                       54698, 37394, 54782,
                                                                       25385, 25475, 38402,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 55874, 0, 3,
                                                                       54782, 37457, 54866,
                                                                       25475, 25565, 38528,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 56042, 0, 3,
                                                                       54866, 37520, 54950,
                                                                       25565, 25655, 38654,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 56210, 0, 3,
                                                                       55034, 37898, 55202,
                                                                       25835, 25985, 39200,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 56490, 0, 3,
                                                                       55202, 38024, 55370,
                                                                       25985, 26135, 39410,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 56770, 0, 3,
                                                                       55370, 38150, 55538,
                                                                       26135, 26285, 39620,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 57050, 0, 3,
                                                                       55538, 38276, 55706,
                                                                       26285, 26435, 39830,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 57330, 0, 3,
                                                                       55706, 38402, 55874,
                                                                       26435, 26585, 40040,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 57610, 0, 3,
                                                                       55874, 38528, 56042,
                                                                       26585, 26735, 40250,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 57890, 0, 3,
                                                                       56210, 39200, 56490,
                                                                       27035, 27260, 41090,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 58310, 0, 3,
                                                                       56490, 39410, 56770,
                                                                       27260, 27485, 41405,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 58730, 0, 3,
                                                                       56770, 39620, 57050,
                                                                       27485, 27710, 41720,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 59150, 0, 3,
                                                                       57050, 39830, 57330,
                                                                       27710, 27935, 42035,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 59570, 0, 3,
                                                                       57330, 40040, 57610,
                                                                       27935, 28160, 42350,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 59990, 0, 3,
                                                                       57890, 41090, 58310,
                                                                       28610, 28925, 43547,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 60578, 0, 3,
                                                                       58310, 41405, 58730,
                                                                       28925, 29240, 43988,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 61166, 0, 3,
                                                                       58730, 41720, 59150,
                                                                       29240, 29555, 44429,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 61754, 0, 3,
                                                                       59150, 42035, 59570,
                                                                       29555, 29870, 44870,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 62342, 0, 3,
                                                                       59990, 43547, 60578,
                                                                       30500, 30920, 46487,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 63126, 0, 3,
                                                                       60578, 43988, 61166,
                                                                       30920, 31340, 47075,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 63910, 0, 3,
                                                                       61166, 44429, 61754,
                                                                       31340, 31760, 47663,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 64694, 0, 3,
                                                                       62342, 46487, 63126,
                                                                       32600, 33140, 49763,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 65702, 0, 3,
                                                                       63126, 47075, 63910,
                                                                       33140, 33680, 50519,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 66710, 0, 3,
                                                                       64694, 49763, 65702,
                                                                       34760, 35435, 53165,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 67970, 3, 36785,
                                                                       36806, 54110, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 68006, 3, 36806,
                                                                       36827, 54138, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 68042, 3, 36827,
                                                                       36848, 54166, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 68078, 3, 36848,
                                                                       36869, 54194, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 68114, 3, 36869,
                                                                       36890, 54222, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 68150, 3, 36890,
                                                                       36911, 54250, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 68186, 3, 36911,
                                                                       36932, 54278, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 68222, 3, 36932,
                                                                       36953, 54306, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 68258, 3, 36953,
                                                                       36974, 54334, ncols,
                                                                       gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 68294, 0, 3,
                                                                       67970, 54110, 68006,
                                                                       37016, 37079, 54362,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 68402, 0, 3,
                                                                       68006, 54138, 68042,
                                                                       37079, 37142, 54446,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 68510, 0, 3,
                                                                       68042, 54166, 68078,
                                                                       37142, 37205, 54530,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 68618, 0, 3,
                                                                       68078, 54194, 68114,
                                                                       37205, 37268, 54614,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 68726, 0, 3,
                                                                       68114, 54222, 68150,
                                                                       37268, 37331, 54698,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 68834, 0, 3,
                                                                       68150, 54250, 68186,
                                                                       37331, 37394, 54782,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 68942, 0, 3,
                                                                       68186, 54278, 68222,
                                                                       37394, 37457, 54866,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 69050, 0, 3,
                                                                       68222, 54306, 68258,
                                                                       37457, 37520, 54950,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 69158, 0, 3,
                                                                       68294, 54362, 68402,
                                                                       37646, 37772, 55034,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 69374, 0, 3,
                                                                       68402, 54446, 68510,
                                                                       37772, 37898, 55202,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 69590, 0, 3,
                                                                       68510, 54530, 68618,
                                                                       37898, 38024, 55370,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 69806, 0, 3,
                                                                       68618, 54614, 68726,
                                                                       38024, 38150, 55538,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 70022, 0, 3,
                                                                       68726, 54698, 68834,
                                                                       38150, 38276, 55706,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 70238, 0, 3,
                                                                       68834, 54782, 68942,
                                                                       38276, 38402, 55874,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 70454, 0, 3,
                                                                       68942, 54866, 69050,
                                                                       38402, 38528, 56042,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 70670, 0, 3,
                                                                       69158, 55034, 69374,
                                                                       38780, 38990, 56210,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 71030, 0, 3,
                                                                       69374, 55202, 69590,
                                                                       38990, 39200, 56490,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 71390, 0, 3,
                                                                       69590, 55370, 69806,
                                                                       39200, 39410, 56770,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 71750, 0, 3,
                                                                       69806, 55538, 70022,
                                                                       39410, 39620, 57050,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 72110, 0, 3,
                                                                       70022, 55706, 70238,
                                                                       39620, 39830, 57330,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 72470, 0, 3,
                                                                       70238, 55874, 70454,
                                                                       39830, 40040, 57610,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 72830, 0, 3,
                                                                       70670, 56210, 71030,
                                                                       40460, 40775, 57890,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 73370, 0, 3,
                                                                       71030, 56490, 71390,
                                                                       40775, 41090, 58310,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 73910, 0, 3,
                                                                       71390, 56770, 71750,
                                                                       41090, 41405, 58730,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 74450, 0, 3,
                                                                       71750, 57050, 72110,
                                                                       41405, 41720, 59150,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 74990, 0, 3,
                                                                       72110, 57330, 72470,
                                                                       41720, 42035, 59570,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 75530, 0, 3,
                                                                       72830, 57890, 73370,
                                                                       42665, 43106, 59990,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 76286, 0, 3,
                                                                       73370, 58310, 73910,
                                                                       43106, 43547, 60578,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 77042, 0, 3,
                                                                       73910, 58730, 74450,
                                                                       43547, 43988, 61166,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 77798, 0, 3,
                                                                       74450, 59150, 74990,
                                                                       43988, 44429, 61754,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 78554, 0, 3,
                                                                       75530, 59990, 76286,
                                                                       45311, 45899, 62342,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 79562, 0, 3,
                                                                       76286, 60578, 77042,
                                                                       45899, 46487, 63126,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 80570, 0, 3,
                                                                       77042, 61166, 77798,
                                                                       46487, 47075, 63910,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 81578, 0, 3,
                                                                       78554, 62342, 79562,
                                                                       48251, 49007, 64694,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 82874, 0, 3,
                                                                       79562, 63126, 80570,
                                                                       49007, 49763, 65702,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 84170, 0, 3,
                                                                       81578, 64694, 82874,
                                                                       51275, 52220, 66710,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 85790, 78554, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 87218, 81578, 1296, ncols);

                    simdfunc::contract_primitives(buffer, 89054, 84170, 1620, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 86798, 85790, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 88514, 87218, 36, 1, nmax);

        simdtrf::transform_k_inner(buffer, 90674, 89054, 45, 1, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 91349, 86798, 88514, 15, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 92609, 88514, 90674, 15, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 94229, 91349, 92609, 15, nmax);

        simdtrf::transform_i_inner(buffer, 96749, 94229, 6, 15, nmax);

        simdtrf::transform_d_outer(values + n * npairs, nvalues, buffer, 96749, 195, nmax);
    }

    for (size_t m = 0; m < 975; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
