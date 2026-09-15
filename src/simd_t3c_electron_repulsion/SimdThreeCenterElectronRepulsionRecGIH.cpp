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


#include "SimdThreeCenterElectronRepulsionRecGIH.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSHD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSID.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSND.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNS.hpp"
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
#include "SimdTransferDI.hpp"
#include "SimdTransferDK.hpp"
#include "SimdTransferDL.hpp"
#include "SimdTransferFI.hpp"
#include "SimdTransferFK.hpp"
#include "SimdTransferGI.hpp"
#include "SimdTransferPI.hpp"
#include "SimdTransferPK.hpp"
#include "SimdTransferPL.hpp"
#include "SimdTransferPM.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_gih_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_gih_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 103545, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1287 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 103545, 69774, 6634, dimensions);

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

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1493, 0, 3, 890,
                                                                       926, 1178, 1223, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1548, 0, 3, 926,
                                                                       962, 1223, 1268, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1603, 0, 3, 962,
                                                                       998, 1268, 1313, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1658, 0, 3, 998,
                                                                       1034, 1313, 1358, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1713, 0, 3, 1034,
                                                                       1070, 1358, 1403, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1768, 0, 3, 1070,
                                                                       1106, 1403, 1448, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 1823, 0, 3, 1178,
                                                                       1223, 1493, 1548, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 1889, 0, 3, 1223,
                                                                       1268, 1548, 1603, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 1955, 0, 3, 1268,
                                                                       1313, 1603, 1658, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2021, 0, 3, 1313,
                                                                       1358, 1658, 1713, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2087, 0, 3, 1358,
                                                                       1403, 1713, 1768, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2153, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2156, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2159, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2162, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2165, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2168, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2171, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2174, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2177, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2180, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2183, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2186, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2189, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2192, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2195, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2198, 3, 10, 29,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2207, 3, 11, 32,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2216, 3, 12, 35,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2225, 3, 13, 38,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2234, 3, 14, 41,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2243, 3, 15, 44,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2252, 3, 16, 47,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2261, 3, 17, 50,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2270, 3, 18, 53,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2279, 3, 19, 56,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2288, 3, 20, 59,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2297, 3, 21, 62,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2306, 3, 23, 65,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2324, 3, 26, 71,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2342, 3, 29, 77,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2360, 3, 32, 83,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2378, 3, 35, 89,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2396, 3, 38, 95,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2414, 3, 41, 101,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2432, 3, 44, 107,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2450, 3, 47, 113,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2468, 3, 50, 119,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2486, 3, 53, 125,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2504, 3, 56, 131,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2522, 3, 59, 137,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2540, 3, 65, 143,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2570, 3, 71, 153,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2600, 3, 77, 163,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2630, 3, 83, 173,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2660, 3, 89, 183,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2690, 3, 95, 193,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2720, 3, 101, 203,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2750, 3, 107, 213,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2780, 3, 113, 223,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2810, 3, 119, 233,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2840, 3, 125, 243,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2870, 3, 131, 253,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2900, 3, 143, 263,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2945, 3, 153, 278,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2990, 3, 163, 293,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3035, 3, 173, 308,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3080, 3, 183, 323,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3125, 3, 193, 338,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3170, 3, 203, 353,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3215, 3, 213, 368,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3260, 3, 223, 383,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3305, 3, 233, 398,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3350, 3, 243, 413,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3395, 3, 263, 428,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3458, 3, 278, 449,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3521, 3, 293, 470,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3584, 3, 308, 491,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3647, 3, 323, 512,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3710, 3, 338, 533,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3773, 3, 353, 554,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3836, 3, 368, 575,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3899, 3, 383, 596,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3962, 3, 398, 617,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4025, 3, 428, 638,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4109, 3, 449, 666,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4193, 3, 470, 694,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4277, 3, 491, 722,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4361, 3, 512, 750,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4445, 3, 533, 778,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4529, 3, 554, 806,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4613, 3, 575, 834,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4697, 3, 596, 862,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4781, 3, 638, 890,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4889, 3, 666, 926,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4997, 3, 694, 962,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5105, 3, 722, 998,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5213, 3, 750,
                                                                       1034, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5321, 3, 778,
                                                                       1070, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5429, 3, 806,
                                                                       1106, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5537, 3, 834,
                                                                       1142, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5645, 3, 890,
                                                                       1178, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5780, 3, 926,
                                                                       1223, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5915, 3, 962,
                                                                       1268, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6050, 3, 998,
                                                                       1313, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6185, 3, 1034,
                                                                       1358, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6320, 3, 1070,
                                                                       1403, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6455, 3, 1106,
                                                                       1448, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 6590, 3, 1178,
                                                                       1493, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 6755, 3, 1223,
                                                                       1548, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 6920, 3, 1268,
                                                                       1603, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 7085, 3, 1313,
                                                                       1658, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 7250, 3, 1358,
                                                                       1713, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 7415, 3, 1403,
                                                                       1768, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 7580, 3, 1493,
                                                                       1823, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 7778, 3, 1548,
                                                                       1889, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 7976, 3, 1603,
                                                                       1955, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 8174, 3, 1658,
                                                                       2021, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 8372, 3, 1713,
                                                                       2087, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8570, 3, 8, 9,
                                                                       2159, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8576, 3, 9, 10,
                                                                       2162, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8582, 3, 10, 11,
                                                                       2165, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8588, 3, 11, 12,
                                                                       2168, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8594, 3, 12, 13,
                                                                       2171, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8600, 3, 13, 14,
                                                                       2174, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8606, 3, 14, 15,
                                                                       2177, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8612, 3, 15, 16,
                                                                       2180, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8618, 3, 16, 17,
                                                                       2183, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8624, 3, 17, 18,
                                                                       2186, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8630, 3, 18, 19,
                                                                       2189, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8636, 3, 19, 20,
                                                                       2192, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8642, 3, 20, 21,
                                                                       2195, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8648, 0, 3, 8570,
                                                                       2159, 8576, 2198, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8666, 0, 3, 8576,
                                                                       2162, 8582, 2207, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8684, 0, 3, 8582,
                                                                       2165, 8588, 2216, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8702, 0, 3, 8588,
                                                                       2168, 8594, 2225, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8720, 0, 3, 8594,
                                                                       2171, 8600, 2234, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8738, 0, 3, 8600,
                                                                       2174, 8606, 2243, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8756, 0, 3, 8606,
                                                                       2177, 8612, 2252, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8774, 0, 3, 8612,
                                                                       2180, 8618, 2261, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8792, 0, 3, 8618,
                                                                       2183, 8624, 2270, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8810, 0, 3, 8624,
                                                                       2186, 8630, 2279, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8828, 0, 3, 8630,
                                                                       2189, 8636, 2288, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8846, 0, 3, 8636,
                                                                       2192, 8642, 2297, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8864, 0, 3, 8648,
                                                                       2198, 8666, 65, 71, 2342,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8900, 0, 3, 8666,
                                                                       2207, 8684, 71, 77, 2360,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8936, 0, 3, 8684,
                                                                       2216, 8702, 77, 83, 2378,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8972, 0, 3, 8702,
                                                                       2225, 8720, 83, 89, 2396,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9008, 0, 3, 8720,
                                                                       2234, 8738, 89, 95, 2414,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9044, 0, 3, 8738,
                                                                       2243, 8756, 95, 101, 2432,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9080, 0, 3, 8756,
                                                                       2252, 8774, 101, 107,
                                                                       2450, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9116, 0, 3, 8774,
                                                                       2261, 8792, 107, 113,
                                                                       2468, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9152, 0, 3, 8792,
                                                                       2270, 8810, 113, 119,
                                                                       2486, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9188, 0, 3, 8810,
                                                                       2279, 8828, 119, 125,
                                                                       2504, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9224, 0, 3, 8828,
                                                                       2288, 8846, 125, 131,
                                                                       2522, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9260, 0, 3, 8864,
                                                                       2342, 8900, 143, 153,
                                                                       2600, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9320, 0, 3, 8900,
                                                                       2360, 8936, 153, 163,
                                                                       2630, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9380, 0, 3, 8936,
                                                                       2378, 8972, 163, 173,
                                                                       2660, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9440, 0, 3, 8972,
                                                                       2396, 9008, 173, 183,
                                                                       2690, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9500, 0, 3, 9008,
                                                                       2414, 9044, 183, 193,
                                                                       2720, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9560, 0, 3, 9044,
                                                                       2432, 9080, 193, 203,
                                                                       2750, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9620, 0, 3, 9080,
                                                                       2450, 9116, 203, 213,
                                                                       2780, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9680, 0, 3, 9116,
                                                                       2468, 9152, 213, 223,
                                                                       2810, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9740, 0, 3, 9152,
                                                                       2486, 9188, 223, 233,
                                                                       2840, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9800, 0, 3, 9188,
                                                                       2504, 9224, 233, 243,
                                                                       2870, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9860, 0, 3, 9260,
                                                                       2600, 9320, 263, 278,
                                                                       2990, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9950, 0, 3, 9320,
                                                                       2630, 9380, 278, 293,
                                                                       3035, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10040, 0, 3, 9380,
                                                                       2660, 9440, 293, 308,
                                                                       3080, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10130, 0, 3, 9440,
                                                                       2690, 9500, 308, 323,
                                                                       3125, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10220, 0, 3, 9500,
                                                                       2720, 9560, 323, 338,
                                                                       3170, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10310, 0, 3, 9560,
                                                                       2750, 9620, 338, 353,
                                                                       3215, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10400, 0, 3, 9620,
                                                                       2780, 9680, 353, 368,
                                                                       3260, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10490, 0, 3, 9680,
                                                                       2810, 9740, 368, 383,
                                                                       3305, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10580, 0, 3, 9740,
                                                                       2840, 9800, 383, 398,
                                                                       3350, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10670, 0, 3, 9860,
                                                                       2990, 9950, 428, 449,
                                                                       3521, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10796, 0, 3, 9950,
                                                                       3035, 10040, 449, 470,
                                                                       3584, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10922, 0, 3,
                                                                       10040, 3080, 10130, 470,
                                                                       491, 3647, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11048, 0, 3,
                                                                       10130, 3125, 10220, 491,
                                                                       512, 3710, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11174, 0, 3,
                                                                       10220, 3170, 10310, 512,
                                                                       533, 3773, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11300, 0, 3,
                                                                       10310, 3215, 10400, 533,
                                                                       554, 3836, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11426, 0, 3,
                                                                       10400, 3260, 10490, 554,
                                                                       575, 3899, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11552, 0, 3,
                                                                       10490, 3305, 10580, 575,
                                                                       596, 3962, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11678, 0, 3,
                                                                       10670, 3521, 10796, 638,
                                                                       666, 4193, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11846, 0, 3,
                                                                       10796, 3584, 10922, 666,
                                                                       694, 4277, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 12014, 0, 3,
                                                                       10922, 3647, 11048, 694,
                                                                       722, 4361, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 12182, 0, 3,
                                                                       11048, 3710, 11174, 722,
                                                                       750, 4445, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 12350, 0, 3,
                                                                       11174, 3773, 11300, 750,
                                                                       778, 4529, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 12518, 0, 3,
                                                                       11300, 3836, 11426, 778,
                                                                       806, 4613, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 12686, 0, 3,
                                                                       11426, 3899, 11552, 806,
                                                                       834, 4697, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 12854, 0, 3,
                                                                       11678, 4193, 11846, 890,
                                                                       926, 4997, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 13070, 0, 3,
                                                                       11846, 4277, 12014, 926,
                                                                       962, 5105, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 13286, 0, 3,
                                                                       12014, 4361, 12182, 962,
                                                                       998, 5213, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 13502, 0, 3,
                                                                       12182, 4445, 12350, 998,
                                                                       1034, 5321, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 13718, 0, 3,
                                                                       12350, 4529, 12518, 1034,
                                                                       1070, 5429, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 13934, 0, 3,
                                                                       12518, 4613, 12686, 1070,
                                                                       1106, 5537, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 14150, 0, 3,
                                                                       12854, 4997, 13070, 1178,
                                                                       1223, 5915, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 14420, 0, 3,
                                                                       13070, 5105, 13286, 1223,
                                                                       1268, 6050, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 14690, 0, 3,
                                                                       13286, 5213, 13502, 1268,
                                                                       1313, 6185, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 14960, 0, 3,
                                                                       13502, 5321, 13718, 1313,
                                                                       1358, 6320, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 15230, 0, 3,
                                                                       13718, 5429, 13934, 1358,
                                                                       1403, 6455, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 15500, 0, 3,
                                                                       14150, 5915, 14420, 1493,
                                                                       1548, 6920, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 15830, 0, 3,
                                                                       14420, 6050, 14690, 1548,
                                                                       1603, 7085, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 16160, 0, 3,
                                                                       14690, 6185, 14960, 1603,
                                                                       1658, 7250, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 16490, 0, 3,
                                                                       14960, 6320, 15230, 1658,
                                                                       1713, 7415, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 16820, 0, 3,
                                                                       15500, 6920, 15830, 1823,
                                                                       1889, 7976, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 17216, 0, 3,
                                                                       15830, 7085, 16160, 1889,
                                                                       1955, 8174, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 17612, 0, 3,
                                                                       16160, 7250, 16490, 1955,
                                                                       2021, 8372, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18008, 3, 2153,
                                                                       2156, 8570, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18018, 3, 2156,
                                                                       2159, 8576, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18028, 3, 2159,
                                                                       2162, 8582, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18038, 3, 2162,
                                                                       2165, 8588, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18048, 3, 2165,
                                                                       2168, 8594, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18058, 3, 2168,
                                                                       2171, 8600, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18068, 3, 2171,
                                                                       2174, 8606, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18078, 3, 2174,
                                                                       2177, 8612, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18088, 3, 2177,
                                                                       2180, 8618, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18098, 3, 2180,
                                                                       2183, 8624, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18108, 3, 2183,
                                                                       2186, 8630, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18118, 3, 2186,
                                                                       2189, 8636, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18128, 3, 2189,
                                                                       2192, 8642, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18138, 0, 3,
                                                                       18008, 8570, 18018, 8648,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18168, 0, 3,
                                                                       18018, 8576, 18028, 8666,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18198, 0, 3,
                                                                       18028, 8582, 18038, 8684,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18228, 0, 3,
                                                                       18038, 8588, 18048, 8702,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18258, 0, 3,
                                                                       18048, 8594, 18058, 8720,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18288, 0, 3,
                                                                       18058, 8600, 18068, 8738,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18318, 0, 3,
                                                                       18068, 8606, 18078, 8756,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18348, 0, 3,
                                                                       18078, 8612, 18088, 8774,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18378, 0, 3,
                                                                       18088, 8618, 18098, 8792,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18408, 0, 3,
                                                                       18098, 8624, 18108, 8810,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18438, 0, 3,
                                                                       18108, 8630, 18118, 8828,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18468, 0, 3,
                                                                       18118, 8636, 18128, 8846,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18498, 0, 3,
                                                                       18138, 8648, 18168, 2306,
                                                                       2324, 8864, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18558, 0, 3,
                                                                       18168, 8666, 18198, 2324,
                                                                       2342, 8900, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18618, 0, 3,
                                                                       18198, 8684, 18228, 2342,
                                                                       2360, 8936, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18678, 0, 3,
                                                                       18228, 8702, 18258, 2360,
                                                                       2378, 8972, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18738, 0, 3,
                                                                       18258, 8720, 18288, 2378,
                                                                       2396, 9008, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18798, 0, 3,
                                                                       18288, 8738, 18318, 2396,
                                                                       2414, 9044, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18858, 0, 3,
                                                                       18318, 8756, 18348, 2414,
                                                                       2432, 9080, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18918, 0, 3,
                                                                       18348, 8774, 18378, 2432,
                                                                       2450, 9116, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18978, 0, 3,
                                                                       18378, 8792, 18408, 2450,
                                                                       2468, 9152, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 19038, 0, 3,
                                                                       18408, 8810, 18438, 2468,
                                                                       2486, 9188, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 19098, 0, 3,
                                                                       18438, 8828, 18468, 2486,
                                                                       2504, 9224, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19158, 0, 3,
                                                                       18498, 8864, 18558, 2540,
                                                                       2570, 9260, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19258, 0, 3,
                                                                       18558, 8900, 18618, 2570,
                                                                       2600, 9320, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19358, 0, 3,
                                                                       18618, 8936, 18678, 2600,
                                                                       2630, 9380, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19458, 0, 3,
                                                                       18678, 8972, 18738, 2630,
                                                                       2660, 9440, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19558, 0, 3,
                                                                       18738, 9008, 18798, 2660,
                                                                       2690, 9500, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19658, 0, 3,
                                                                       18798, 9044, 18858, 2690,
                                                                       2720, 9560, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19758, 0, 3,
                                                                       18858, 9080, 18918, 2720,
                                                                       2750, 9620, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19858, 0, 3,
                                                                       18918, 9116, 18978, 2750,
                                                                       2780, 9680, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19958, 0, 3,
                                                                       18978, 9152, 19038, 2780,
                                                                       2810, 9740, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 20058, 0, 3,
                                                                       19038, 9188, 19098, 2810,
                                                                       2840, 9800, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 20158, 0, 3,
                                                                       19158, 9260, 19258, 2900,
                                                                       2945, 9860, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 20308, 0, 3,
                                                                       19258, 9320, 19358, 2945,
                                                                       2990, 9950, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 20458, 0, 3,
                                                                       19358, 9380, 19458, 2990,
                                                                       3035, 10040, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 20608, 0, 3,
                                                                       19458, 9440, 19558, 3035,
                                                                       3080, 10130, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 20758, 0, 3,
                                                                       19558, 9500, 19658, 3080,
                                                                       3125, 10220, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 20908, 0, 3,
                                                                       19658, 9560, 19758, 3125,
                                                                       3170, 10310, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 21058, 0, 3,
                                                                       19758, 9620, 19858, 3170,
                                                                       3215, 10400, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 21208, 0, 3,
                                                                       19858, 9680, 19958, 3215,
                                                                       3260, 10490, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 21358, 0, 3,
                                                                       19958, 9740, 20058, 3260,
                                                                       3305, 10580, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 21508, 0, 3,
                                                                       20158, 9860, 20308, 3395,
                                                                       3458, 10670, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 21718, 0, 3,
                                                                       20308, 9950, 20458, 3458,
                                                                       3521, 10796, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 21928, 0, 3,
                                                                       20458, 10040, 20608, 3521,
                                                                       3584, 10922, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 22138, 0, 3,
                                                                       20608, 10130, 20758, 3584,
                                                                       3647, 11048, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 22348, 0, 3,
                                                                       20758, 10220, 20908, 3647,
                                                                       3710, 11174, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 22558, 0, 3,
                                                                       20908, 10310, 21058, 3710,
                                                                       3773, 11300, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 22768, 0, 3,
                                                                       21058, 10400, 21208, 3773,
                                                                       3836, 11426, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 22978, 0, 3,
                                                                       21208, 10490, 21358, 3836,
                                                                       3899, 11552, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 23188, 0, 3,
                                                                       21508, 10670, 21718, 4025,
                                                                       4109, 11678, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 23468, 0, 3,
                                                                       21718, 10796, 21928, 4109,
                                                                       4193, 11846, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 23748, 0, 3,
                                                                       21928, 10922, 22138, 4193,
                                                                       4277, 12014, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 24028, 0, 3,
                                                                       22138, 11048, 22348, 4277,
                                                                       4361, 12182, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 24308, 0, 3,
                                                                       22348, 11174, 22558, 4361,
                                                                       4445, 12350, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 24588, 0, 3,
                                                                       22558, 11300, 22768, 4445,
                                                                       4529, 12518, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 24868, 0, 3,
                                                                       22768, 11426, 22978, 4529,
                                                                       4613, 12686, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 25148, 0, 3,
                                                                       23188, 11678, 23468, 4781,
                                                                       4889, 12854, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 25508, 0, 3,
                                                                       23468, 11846, 23748, 4889,
                                                                       4997, 13070, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 25868, 0, 3,
                                                                       23748, 12014, 24028, 4997,
                                                                       5105, 13286, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 26228, 0, 3,
                                                                       24028, 12182, 24308, 5105,
                                                                       5213, 13502, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 26588, 0, 3,
                                                                       24308, 12350, 24588, 5213,
                                                                       5321, 13718, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 26948, 0, 3,
                                                                       24588, 12518, 24868, 5321,
                                                                       5429, 13934, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 27308, 0, 3,
                                                                       25148, 12854, 25508, 5645,
                                                                       5780, 14150, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 27758, 0, 3,
                                                                       25508, 13070, 25868, 5780,
                                                                       5915, 14420, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 28208, 0, 3,
                                                                       25868, 13286, 26228, 5915,
                                                                       6050, 14690, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 28658, 0, 3,
                                                                       26228, 13502, 26588, 6050,
                                                                       6185, 14960, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 29108, 0, 3,
                                                                       26588, 13718, 26948, 6185,
                                                                       6320, 15230, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 29558, 0, 3,
                                                                       27308, 14150, 27758, 6590,
                                                                       6755, 15500, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 30108, 0, 3,
                                                                       27758, 14420, 28208, 6755,
                                                                       6920, 15830, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 30658, 0, 3,
                                                                       28208, 14690, 28658, 6920,
                                                                       7085, 16160, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 31208, 0, 3,
                                                                       28658, 14960, 29108, 7085,
                                                                       7250, 16490, ncols, gamma,
                                                                       p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 31758, 0, 3,
                                                                       29558, 15500, 30108, 7580,
                                                                       7778, 16820, ncols, gamma,
                                                                       p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 32418, 0, 3,
                                                                       30108, 15830, 30658, 7778,
                                                                       7976, 17216, ncols, gamma,
                                                                       p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 33078, 0, 3,
                                                                       30658, 16160, 31208, 7976,
                                                                       8174, 17612, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33738, 3, 8570,
                                                                       8576, 18028, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33753, 3, 8576,
                                                                       8582, 18038, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33768, 3, 8582,
                                                                       8588, 18048, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33783, 3, 8588,
                                                                       8594, 18058, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33798, 3, 8594,
                                                                       8600, 18068, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33813, 3, 8600,
                                                                       8606, 18078, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33828, 3, 8606,
                                                                       8612, 18088, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33843, 3, 8612,
                                                                       8618, 18098, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33858, 3, 8618,
                                                                       8624, 18108, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33873, 3, 8624,
                                                                       8630, 18118, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33888, 3, 8630,
                                                                       8636, 18128, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 33903, 0, 3,
                                                                       33738, 18028, 33753, 8648,
                                                                       8666, 18198, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 33948, 0, 3,
                                                                       33753, 18038, 33768, 8666,
                                                                       8684, 18228, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 33993, 0, 3,
                                                                       33768, 18048, 33783, 8684,
                                                                       8702, 18258, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 34038, 0, 3,
                                                                       33783, 18058, 33798, 8702,
                                                                       8720, 18288, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 34083, 0, 3,
                                                                       33798, 18068, 33813, 8720,
                                                                       8738, 18318, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 34128, 0, 3,
                                                                       33813, 18078, 33828, 8738,
                                                                       8756, 18348, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 34173, 0, 3,
                                                                       33828, 18088, 33843, 8756,
                                                                       8774, 18378, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 34218, 0, 3,
                                                                       33843, 18098, 33858, 8774,
                                                                       8792, 18408, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 34263, 0, 3,
                                                                       33858, 18108, 33873, 8792,
                                                                       8810, 18438, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 34308, 0, 3,
                                                                       33873, 18118, 33888, 8810,
                                                                       8828, 18468, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 34353, 0, 3,
                                                                       33903, 18198, 33948, 8864,
                                                                       8900, 18618, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 34443, 0, 3,
                                                                       33948, 18228, 33993, 8900,
                                                                       8936, 18678, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 34533, 0, 3,
                                                                       33993, 18258, 34038, 8936,
                                                                       8972, 18738, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 34623, 0, 3,
                                                                       34038, 18288, 34083, 8972,
                                                                       9008, 18798, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 34713, 0, 3,
                                                                       34083, 18318, 34128, 9008,
                                                                       9044, 18858, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 34803, 0, 3,
                                                                       34128, 18348, 34173, 9044,
                                                                       9080, 18918, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 34893, 0, 3,
                                                                       34173, 18378, 34218, 9080,
                                                                       9116, 18978, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 34983, 0, 3,
                                                                       34218, 18408, 34263, 9116,
                                                                       9152, 19038, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 35073, 0, 3,
                                                                       34263, 18438, 34308, 9152,
                                                                       9188, 19098, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 35163, 0, 3,
                                                                       34353, 18618, 34443, 9260,
                                                                       9320, 19358, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 35313, 0, 3,
                                                                       34443, 18678, 34533, 9320,
                                                                       9380, 19458, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 35463, 0, 3,
                                                                       34533, 18738, 34623, 9380,
                                                                       9440, 19558, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 35613, 0, 3,
                                                                       34623, 18798, 34713, 9440,
                                                                       9500, 19658, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 35763, 0, 3,
                                                                       34713, 18858, 34803, 9500,
                                                                       9560, 19758, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 35913, 0, 3,
                                                                       34803, 18918, 34893, 9560,
                                                                       9620, 19858, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 36063, 0, 3,
                                                                       34893, 18978, 34983, 9620,
                                                                       9680, 19958, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 36213, 0, 3,
                                                                       34983, 19038, 35073, 9680,
                                                                       9740, 20058, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 36363, 0, 3,
                                                                       35163, 19358, 35313, 9860,
                                                                       9950, 20458, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 36588, 0, 3,
                                                                       35313, 19458, 35463, 9950,
                                                                       10040, 20608, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 36813, 0, 3,
                                                                       35463, 19558, 35613,
                                                                       10040, 10130, 20758,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 37038, 0, 3,
                                                                       35613, 19658, 35763,
                                                                       10130, 10220, 20908,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 37263, 0, 3,
                                                                       35763, 19758, 35913,
                                                                       10220, 10310, 21058,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 37488, 0, 3,
                                                                       35913, 19858, 36063,
                                                                       10310, 10400, 21208,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 37713, 0, 3,
                                                                       36063, 19958, 36213,
                                                                       10400, 10490, 21358,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 37938, 0, 3,
                                                                       36363, 20458, 36588,
                                                                       10670, 10796, 21928,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 38253, 0, 3,
                                                                       36588, 20608, 36813,
                                                                       10796, 10922, 22138,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 38568, 0, 3,
                                                                       36813, 20758, 37038,
                                                                       10922, 11048, 22348,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 38883, 0, 3,
                                                                       37038, 20908, 37263,
                                                                       11048, 11174, 22558,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 39198, 0, 3,
                                                                       37263, 21058, 37488,
                                                                       11174, 11300, 22768,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 39513, 0, 3,
                                                                       37488, 21208, 37713,
                                                                       11300, 11426, 22978,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 39828, 0, 3,
                                                                       37938, 21928, 38253,
                                                                       11678, 11846, 23748,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 40248, 0, 3,
                                                                       38253, 22138, 38568,
                                                                       11846, 12014, 24028,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 40668, 0, 3,
                                                                       38568, 22348, 38883,
                                                                       12014, 12182, 24308,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 41088, 0, 3,
                                                                       38883, 22558, 39198,
                                                                       12182, 12350, 24588,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 41508, 0, 3,
                                                                       39198, 22768, 39513,
                                                                       12350, 12518, 24868,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 41928, 0, 3,
                                                                       39828, 23748, 40248,
                                                                       12854, 13070, 25868,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 42468, 0, 3,
                                                                       40248, 24028, 40668,
                                                                       13070, 13286, 26228,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 43008, 0, 3,
                                                                       40668, 24308, 41088,
                                                                       13286, 13502, 26588,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 43548, 0, 3,
                                                                       41088, 24588, 41508,
                                                                       13502, 13718, 26948,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 44088, 0, 3,
                                                                       41928, 25868, 42468,
                                                                       14150, 14420, 28208,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 44763, 0, 3,
                                                                       42468, 26228, 43008,
                                                                       14420, 14690, 28658,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 45438, 0, 3,
                                                                       43008, 26588, 43548,
                                                                       14690, 14960, 29108,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 46113, 0, 3,
                                                                       44088, 28208, 44763,
                                                                       15500, 15830, 30658,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 46938, 0, 3,
                                                                       44763, 28658, 45438,
                                                                       15830, 16160, 31208,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 47763, 0, 3,
                                                                       46113, 30658, 46938,
                                                                       16820, 17216, 33078,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 48753, 3, 18008,
                                                                       18018, 33738, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 48774, 3, 18018,
                                                                       18028, 33753, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 48795, 3, 18028,
                                                                       18038, 33768, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 48816, 3, 18038,
                                                                       18048, 33783, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 48837, 3, 18048,
                                                                       18058, 33798, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 48858, 3, 18058,
                                                                       18068, 33813, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 48879, 3, 18068,
                                                                       18078, 33828, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 48900, 3, 18078,
                                                                       18088, 33843, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 48921, 3, 18088,
                                                                       18098, 33858, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 48942, 3, 18098,
                                                                       18108, 33873, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 48963, 3, 18108,
                                                                       18118, 33888, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 48984, 0, 3,
                                                                       48753, 33738, 48774,
                                                                       18138, 18168, 33903,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 49047, 0, 3,
                                                                       48774, 33753, 48795,
                                                                       18168, 18198, 33948,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 49110, 0, 3,
                                                                       48795, 33768, 48816,
                                                                       18198, 18228, 33993,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 49173, 0, 3,
                                                                       48816, 33783, 48837,
                                                                       18228, 18258, 34038,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 49236, 0, 3,
                                                                       48837, 33798, 48858,
                                                                       18258, 18288, 34083,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 49299, 0, 3,
                                                                       48858, 33813, 48879,
                                                                       18288, 18318, 34128,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 49362, 0, 3,
                                                                       48879, 33828, 48900,
                                                                       18318, 18348, 34173,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 49425, 0, 3,
                                                                       48900, 33843, 48921,
                                                                       18348, 18378, 34218,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 49488, 0, 3,
                                                                       48921, 33858, 48942,
                                                                       18378, 18408, 34263,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 49551, 0, 3,
                                                                       48942, 33873, 48963,
                                                                       18408, 18438, 34308,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 49614, 0, 3,
                                                                       48984, 33903, 49047,
                                                                       18498, 18558, 34353,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 49740, 0, 3,
                                                                       49047, 33948, 49110,
                                                                       18558, 18618, 34443,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 49866, 0, 3,
                                                                       49110, 33993, 49173,
                                                                       18618, 18678, 34533,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 49992, 0, 3,
                                                                       49173, 34038, 49236,
                                                                       18678, 18738, 34623,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 50118, 0, 3,
                                                                       49236, 34083, 49299,
                                                                       18738, 18798, 34713,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 50244, 0, 3,
                                                                       49299, 34128, 49362,
                                                                       18798, 18858, 34803,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 50370, 0, 3,
                                                                       49362, 34173, 49425,
                                                                       18858, 18918, 34893,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 50496, 0, 3,
                                                                       49425, 34218, 49488,
                                                                       18918, 18978, 34983,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 50622, 0, 3,
                                                                       49488, 34263, 49551,
                                                                       18978, 19038, 35073,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 50748, 0, 3,
                                                                       49614, 34353, 49740,
                                                                       19158, 19258, 35163,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 50958, 0, 3,
                                                                       49740, 34443, 49866,
                                                                       19258, 19358, 35313,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 51168, 0, 3,
                                                                       49866, 34533, 49992,
                                                                       19358, 19458, 35463,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 51378, 0, 3,
                                                                       49992, 34623, 50118,
                                                                       19458, 19558, 35613,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 51588, 0, 3,
                                                                       50118, 34713, 50244,
                                                                       19558, 19658, 35763,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 51798, 0, 3,
                                                                       50244, 34803, 50370,
                                                                       19658, 19758, 35913,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 52008, 0, 3,
                                                                       50370, 34893, 50496,
                                                                       19758, 19858, 36063,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 52218, 0, 3,
                                                                       50496, 34983, 50622,
                                                                       19858, 19958, 36213,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 52428, 0, 3,
                                                                       50748, 35163, 50958,
                                                                       20158, 20308, 36363,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 52743, 0, 3,
                                                                       50958, 35313, 51168,
                                                                       20308, 20458, 36588,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 53058, 0, 3,
                                                                       51168, 35463, 51378,
                                                                       20458, 20608, 36813,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 53373, 0, 3,
                                                                       51378, 35613, 51588,
                                                                       20608, 20758, 37038,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 53688, 0, 3,
                                                                       51588, 35763, 51798,
                                                                       20758, 20908, 37263,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 54003, 0, 3,
                                                                       51798, 35913, 52008,
                                                                       20908, 21058, 37488,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 54318, 0, 3,
                                                                       52008, 36063, 52218,
                                                                       21058, 21208, 37713,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 54633, 0, 3,
                                                                       52428, 36363, 52743,
                                                                       21508, 21718, 37938,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 55074, 0, 3,
                                                                       52743, 36588, 53058,
                                                                       21718, 21928, 38253,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 55515, 0, 3,
                                                                       53058, 36813, 53373,
                                                                       21928, 22138, 38568,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 55956, 0, 3,
                                                                       53373, 37038, 53688,
                                                                       22138, 22348, 38883,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 56397, 0, 3,
                                                                       53688, 37263, 54003,
                                                                       22348, 22558, 39198,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 56838, 0, 3,
                                                                       54003, 37488, 54318,
                                                                       22558, 22768, 39513,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 57279, 0, 3,
                                                                       54633, 37938, 55074,
                                                                       23188, 23468, 39828,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 57867, 0, 3,
                                                                       55074, 38253, 55515,
                                                                       23468, 23748, 40248,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 58455, 0, 3,
                                                                       55515, 38568, 55956,
                                                                       23748, 24028, 40668,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 59043, 0, 3,
                                                                       55956, 38883, 56397,
                                                                       24028, 24308, 41088,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 59631, 0, 3,
                                                                       56397, 39198, 56838,
                                                                       24308, 24588, 41508,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 60219, 0, 3,
                                                                       57279, 39828, 57867,
                                                                       25148, 25508, 41928,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 60975, 0, 3,
                                                                       57867, 40248, 58455,
                                                                       25508, 25868, 42468,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 61731, 0, 3,
                                                                       58455, 40668, 59043,
                                                                       25868, 26228, 43008,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 62487, 0, 3,
                                                                       59043, 41088, 59631,
                                                                       26228, 26588, 43548,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 63243, 0, 3,
                                                                       60219, 41928, 60975,
                                                                       27308, 27758, 44088,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 64188, 0, 3,
                                                                       60975, 42468, 61731,
                                                                       27758, 28208, 44763,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 65133, 0, 3,
                                                                       61731, 43008, 62487,
                                                                       28208, 28658, 45438,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 66078, 0, 3,
                                                                       63243, 44088, 64188,
                                                                       29558, 30108, 46113,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 67233, 0, 3,
                                                                       64188, 44763, 65133,
                                                                       30108, 30658, 46938,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 68388, 0, 3,
                                                                       66078, 46113, 67233,
                                                                       31758, 32418, 47763,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 69774, 57279, 588, ncols);

                    simdfunc::contract_primitives(buffer, 70670, 60219, 756, ncols);

                    simdfunc::contract_primitives(buffer, 71822, 63243, 945, ncols);

                    simdfunc::contract_primitives(buffer, 73262, 66078, 1155, ncols);

                    simdfunc::contract_primitives(buffer, 75022, 68388, 1386, ncols);
                }
            }
        }

        simdtrf::transform_h_inner(buffer, 70362, 69774, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 71426, 70670, 36, 1, nmax);

        simdtrf::transform_h_inner(buffer, 72767, 71822, 45, 1, nmax);

        simdtrf::transform_h_inner(buffer, 74417, 73262, 55, 1, nmax);

        simdtrf::transform_h_inner(buffer, 76408, 75022, 66, 1, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 77134, 70362, 71426, 11, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 78058, 71426, 72767, 11, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 79246, 72767, 74417, 11, nmax);

        simdtrf::compute_hrr_pm(buffer, coordinates, 80731, 74417, 76408, 11, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 82546, 77134, 78058, 11, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 84394, 78058, 79246, 11, nmax);

        simdtrf::compute_hrr_dl(buffer, coordinates, 86770, 79246, 80731, 11, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 89740, 82546, 84394, 11, nmax);

        simdtrf::compute_hrr_fk(buffer, coordinates, 92820, 84394, 86770, 11, nmax);

        simdtrf::compute_hrr_gi(buffer, coordinates, 96780, 89740, 92820, 11, nmax);

        simdtrf::transform_i_inner(buffer, 101400, 96780, 15, 11, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 101400, 143, nmax);
    }

    for (size_t m = 0; m < 1287; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
