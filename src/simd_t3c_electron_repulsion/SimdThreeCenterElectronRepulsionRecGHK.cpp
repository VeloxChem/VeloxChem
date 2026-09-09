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


#include "SimdThreeCenterElectronRepulsionRecGHK.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSMD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMS.hpp"
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
#include "SimdTransferDH.hpp"
#include "SimdTransferDI.hpp"
#include "SimdTransferDK.hpp"
#include "SimdTransferFH.hpp"
#include "SimdTransferFI.hpp"
#include "SimdTransferGH.hpp"
#include "SimdTransferPH.hpp"
#include "SimdTransferPI.hpp"
#include "SimdTransferPK.hpp"
#include "SimdTransferPL.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformK.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_ghk_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_ghk_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 158474, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1485 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 158474, 120989, 8610, dimensions);

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

                for (size_t k = 0; k < nprim_c; k++)
                {
                    const auto ncols = dimensions[(i * nprim_b + j) * nprim_c + k];

                    if (ncols == 0) continue;

                    const auto gamma = c_exps[k];

                    const auto q = p + gamma;

                    const auto fq = p * gamma / q;

                    const auto fj = 2.0 * fovl * c_norms[k] * pi * pi * std::sqrt(pi)
                                    / (p * gamma * std::sqrt(q));

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 6, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                                                        16}, ncols, fj, mu, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 23, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 26, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 29, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 32, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 35, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 38, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 41, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 44, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 47, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 50, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 53, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 56, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 59, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 62, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 65, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 68, 0, 3, 7, 8,
                                                                       23, 26, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 74, 0, 3, 8, 9,
                                                                       26, 29, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 80, 0, 3, 9, 10,
                                                                       29, 32, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 86, 0, 3, 10, 11,
                                                                       32, 35, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 92, 0, 3, 11, 12,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 98, 0, 3, 12, 13,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 104, 0, 3, 13, 14,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 110, 0, 3, 14, 15,
                                                                       44, 47, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 116, 0, 3, 15, 16,
                                                                       47, 50, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 122, 0, 3, 16, 17,
                                                                       50, 53, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 128, 0, 3, 17, 18,
                                                                       53, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 134, 0, 3, 18, 19,
                                                                       56, 59, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 140, 0, 3, 19, 20,
                                                                       59, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 146, 0, 3, 20, 21,
                                                                       62, 65, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 152, 0, 3, 23, 26,
                                                                       68, 74, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 162, 0, 3, 26, 29,
                                                                       74, 80, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 172, 0, 3, 29, 32,
                                                                       80, 86, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 182, 0, 3, 32, 35,
                                                                       86, 92, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 192, 0, 3, 35, 38,
                                                                       92, 98, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 202, 0, 3, 38, 41,
                                                                       98, 104, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 212, 0, 3, 41, 44,
                                                                       104, 110, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 222, 0, 3, 44, 47,
                                                                       110, 116, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 232, 0, 3, 47, 50,
                                                                       116, 122, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 242, 0, 3, 50, 53,
                                                                       122, 128, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 252, 0, 3, 53, 56,
                                                                       128, 134, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 262, 0, 3, 56, 59,
                                                                       134, 140, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 272, 0, 3, 59, 62,
                                                                       140, 146, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 282, 0, 3, 68, 74,
                                                                       152, 162, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 297, 0, 3, 74, 80,
                                                                       162, 172, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 312, 0, 3, 80, 86,
                                                                       172, 182, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 327, 0, 3, 86, 92,
                                                                       182, 192, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 342, 0, 3, 92, 98,
                                                                       192, 202, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 357, 0, 3, 98,
                                                                       104, 202, 212, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 372, 0, 3, 104,
                                                                       110, 212, 222, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 387, 0, 3, 110,
                                                                       116, 222, 232, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 402, 0, 3, 116,
                                                                       122, 232, 242, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 417, 0, 3, 122,
                                                                       128, 242, 252, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 432, 0, 3, 128,
                                                                       134, 252, 262, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 447, 0, 3, 134,
                                                                       140, 262, 272, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 462, 0, 3, 152,
                                                                       162, 282, 297, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 483, 0, 3, 162,
                                                                       172, 297, 312, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 504, 0, 3, 172,
                                                                       182, 312, 327, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 525, 0, 3, 182,
                                                                       192, 327, 342, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 546, 0, 3, 192,
                                                                       202, 342, 357, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 567, 0, 3, 202,
                                                                       212, 357, 372, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 588, 0, 3, 212,
                                                                       222, 372, 387, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 609, 0, 3, 222,
                                                                       232, 387, 402, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 630, 0, 3, 232,
                                                                       242, 402, 417, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 651, 0, 3, 242,
                                                                       252, 417, 432, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 672, 0, 3, 252,
                                                                       262, 432, 447, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 693, 0, 3, 282,
                                                                       297, 462, 483, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 721, 0, 3, 297,
                                                                       312, 483, 504, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 749, 0, 3, 312,
                                                                       327, 504, 525, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 777, 0, 3, 327,
                                                                       342, 525, 546, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 805, 0, 3, 342,
                                                                       357, 546, 567, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 833, 0, 3, 357,
                                                                       372, 567, 588, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 861, 0, 3, 372,
                                                                       387, 588, 609, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 889, 0, 3, 387,
                                                                       402, 609, 630, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 917, 0, 3, 402,
                                                                       417, 630, 651, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 945, 0, 3, 417,
                                                                       432, 651, 672, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 973, 0, 3, 462,
                                                                       483, 693, 721, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1009, 0, 3, 483,
                                                                       504, 721, 749, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1045, 0, 3, 504,
                                                                       525, 749, 777, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1081, 0, 3, 525,
                                                                       546, 777, 805, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1117, 0, 3, 546,
                                                                       567, 805, 833, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1153, 0, 3, 567,
                                                                       588, 833, 861, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1189, 0, 3, 588,
                                                                       609, 861, 889, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1225, 0, 3, 609,
                                                                       630, 889, 917, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1261, 0, 3, 630,
                                                                       651, 917, 945, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1297, 0, 3, 693,
                                                                       721, 973, 1009, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1342, 0, 3, 721,
                                                                       749, 1009, 1045, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1387, 0, 3, 749,
                                                                       777, 1045, 1081, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1432, 0, 3, 777,
                                                                       805, 1081, 1117, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1477, 0, 3, 805,
                                                                       833, 1117, 1153, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1522, 0, 3, 833,
                                                                       861, 1153, 1189, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1567, 0, 3, 861,
                                                                       889, 1189, 1225, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1612, 0, 3, 889,
                                                                       917, 1225, 1261, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1657, 0, 3, 973,
                                                                       1009, 1297, 1342, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1712, 0, 3, 1009,
                                                                       1045, 1342, 1387, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1767, 0, 3, 1045,
                                                                       1081, 1387, 1432, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1822, 0, 3, 1081,
                                                                       1117, 1432, 1477, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1877, 0, 3, 1117,
                                                                       1153, 1477, 1522, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1932, 0, 3, 1153,
                                                                       1189, 1522, 1567, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1987, 0, 3, 1189,
                                                                       1225, 1567, 1612, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2042, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2045, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2048, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2051, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2054, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2057, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2060, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2063, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2066, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2069, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2072, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2075, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2078, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2081, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2084, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2087, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2090, 3, 9, 29,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2099, 3, 10, 32,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2108, 3, 11, 35,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2117, 3, 12, 38,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2126, 3, 13, 41,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2135, 3, 14, 44,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2144, 3, 15, 47,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2153, 3, 16, 50,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2162, 3, 17, 53,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2171, 3, 18, 56,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2180, 3, 19, 59,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2189, 3, 20, 62,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2198, 3, 21, 65,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2207, 3, 23, 68,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2225, 3, 26, 74,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2243, 3, 29, 80,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2261, 3, 32, 86,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2279, 3, 35, 92,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2297, 3, 38, 98,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2315, 3, 41, 104,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2333, 3, 44, 110,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2351, 3, 47, 116,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2369, 3, 50, 122,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2387, 3, 53, 128,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2405, 3, 56, 134,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2423, 3, 59, 140,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2441, 3, 62, 146,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2459, 3, 68, 152,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2489, 3, 74, 162,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2519, 3, 80, 172,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2549, 3, 86, 182,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2579, 3, 92, 192,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2609, 3, 98, 202,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2639, 3, 104, 212,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2669, 3, 110, 222,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2699, 3, 116, 232,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2729, 3, 122, 242,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2759, 3, 128, 252,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2789, 3, 134, 262,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2819, 3, 140, 272,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2849, 3, 152, 282,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2894, 3, 162, 297,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2939, 3, 172, 312,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2984, 3, 182, 327,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3029, 3, 192, 342,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3074, 3, 202, 357,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3119, 3, 212, 372,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3164, 3, 222, 387,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3209, 3, 232, 402,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3254, 3, 242, 417,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3299, 3, 252, 432,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3344, 3, 262, 447,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3389, 3, 282, 462,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3452, 3, 297, 483,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3515, 3, 312, 504,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3578, 3, 327, 525,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3641, 3, 342, 546,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3704, 3, 357, 567,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3767, 3, 372, 588,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3830, 3, 387, 609,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3893, 3, 402, 630,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3956, 3, 417, 651,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4019, 3, 432, 672,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4082, 3, 462, 693,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4166, 3, 483, 721,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4250, 3, 504, 749,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4334, 3, 525, 777,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4418, 3, 546, 805,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4502, 3, 567, 833,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4586, 3, 588, 861,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4670, 3, 609, 889,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4754, 3, 630, 917,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4838, 3, 651, 945,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4922, 3, 693, 973,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5030, 3, 721,
                                                                       1009, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5138, 3, 749,
                                                                       1045, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5246, 3, 777,
                                                                       1081, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5354, 3, 805,
                                                                       1117, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5462, 3, 833,
                                                                       1153, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5570, 3, 861,
                                                                       1189, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5678, 3, 889,
                                                                       1225, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5786, 3, 917,
                                                                       1261, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5894, 3, 973,
                                                                       1297, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6029, 3, 1009,
                                                                       1342, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6164, 3, 1045,
                                                                       1387, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6299, 3, 1081,
                                                                       1432, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6434, 3, 1117,
                                                                       1477, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6569, 3, 1153,
                                                                       1522, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6704, 3, 1189,
                                                                       1567, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6839, 3, 1225,
                                                                       1612, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 6974, 3, 1297,
                                                                       1657, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 7139, 3, 1342,
                                                                       1712, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 7304, 3, 1387,
                                                                       1767, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 7469, 3, 1432,
                                                                       1822, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 7634, 3, 1477,
                                                                       1877, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 7799, 3, 1522,
                                                                       1932, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 7964, 3, 1567,
                                                                       1987, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8129, 3, 7, 8,
                                                                       2048, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8135, 3, 8, 9,
                                                                       2051, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8141, 3, 9, 10,
                                                                       2054, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8147, 3, 10, 11,
                                                                       2057, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8153, 3, 11, 12,
                                                                       2060, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8159, 3, 12, 13,
                                                                       2063, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8165, 3, 13, 14,
                                                                       2066, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8171, 3, 14, 15,
                                                                       2069, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8177, 3, 15, 16,
                                                                       2072, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8183, 3, 16, 17,
                                                                       2075, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8189, 3, 17, 18,
                                                                       2078, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8195, 3, 18, 19,
                                                                       2081, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8201, 3, 19, 20,
                                                                       2084, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8207, 3, 20, 21,
                                                                       2087, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8213, 0, 3, 8129,
                                                                       2048, 8135, 2090, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8231, 0, 3, 8135,
                                                                       2051, 8141, 2099, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8249, 0, 3, 8141,
                                                                       2054, 8147, 2108, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8267, 0, 3, 8147,
                                                                       2057, 8153, 2117, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8285, 0, 3, 8153,
                                                                       2060, 8159, 2126, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8303, 0, 3, 8159,
                                                                       2063, 8165, 2135, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8321, 0, 3, 8165,
                                                                       2066, 8171, 2144, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8339, 0, 3, 8171,
                                                                       2069, 8177, 2153, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8357, 0, 3, 8177,
                                                                       2072, 8183, 2162, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8375, 0, 3, 8183,
                                                                       2075, 8189, 2171, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8393, 0, 3, 8189,
                                                                       2078, 8195, 2180, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8411, 0, 3, 8195,
                                                                       2081, 8201, 2189, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8429, 0, 3, 8201,
                                                                       2084, 8207, 2198, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8447, 0, 3, 8213,
                                                                       2090, 8231, 68, 74, 2243,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8483, 0, 3, 8231,
                                                                       2099, 8249, 74, 80, 2261,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8519, 0, 3, 8249,
                                                                       2108, 8267, 80, 86, 2279,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8555, 0, 3, 8267,
                                                                       2117, 8285, 86, 92, 2297,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8591, 0, 3, 8285,
                                                                       2126, 8303, 92, 98, 2315,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8627, 0, 3, 8303,
                                                                       2135, 8321, 98, 104, 2333,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8663, 0, 3, 8321,
                                                                       2144, 8339, 104, 110,
                                                                       2351, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8699, 0, 3, 8339,
                                                                       2153, 8357, 110, 116,
                                                                       2369, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8735, 0, 3, 8357,
                                                                       2162, 8375, 116, 122,
                                                                       2387, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8771, 0, 3, 8375,
                                                                       2171, 8393, 122, 128,
                                                                       2405, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8807, 0, 3, 8393,
                                                                       2180, 8411, 128, 134,
                                                                       2423, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8843, 0, 3, 8411,
                                                                       2189, 8429, 134, 140,
                                                                       2441, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8879, 0, 3, 8447,
                                                                       2243, 8483, 152, 162,
                                                                       2519, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8939, 0, 3, 8483,
                                                                       2261, 8519, 162, 172,
                                                                       2549, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8999, 0, 3, 8519,
                                                                       2279, 8555, 172, 182,
                                                                       2579, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9059, 0, 3, 8555,
                                                                       2297, 8591, 182, 192,
                                                                       2609, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9119, 0, 3, 8591,
                                                                       2315, 8627, 192, 202,
                                                                       2639, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9179, 0, 3, 8627,
                                                                       2333, 8663, 202, 212,
                                                                       2669, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9239, 0, 3, 8663,
                                                                       2351, 8699, 212, 222,
                                                                       2699, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9299, 0, 3, 8699,
                                                                       2369, 8735, 222, 232,
                                                                       2729, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9359, 0, 3, 8735,
                                                                       2387, 8771, 232, 242,
                                                                       2759, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9419, 0, 3, 8771,
                                                                       2405, 8807, 242, 252,
                                                                       2789, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9479, 0, 3, 8807,
                                                                       2423, 8843, 252, 262,
                                                                       2819, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9539, 0, 3, 8879,
                                                                       2519, 8939, 282, 297,
                                                                       2939, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9629, 0, 3, 8939,
                                                                       2549, 8999, 297, 312,
                                                                       2984, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9719, 0, 3, 8999,
                                                                       2579, 9059, 312, 327,
                                                                       3029, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9809, 0, 3, 9059,
                                                                       2609, 9119, 327, 342,
                                                                       3074, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9899, 0, 3, 9119,
                                                                       2639, 9179, 342, 357,
                                                                       3119, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9989, 0, 3, 9179,
                                                                       2669, 9239, 357, 372,
                                                                       3164, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10079, 0, 3, 9239,
                                                                       2699, 9299, 372, 387,
                                                                       3209, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10169, 0, 3, 9299,
                                                                       2729, 9359, 387, 402,
                                                                       3254, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10259, 0, 3, 9359,
                                                                       2759, 9419, 402, 417,
                                                                       3299, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10349, 0, 3, 9419,
                                                                       2789, 9479, 417, 432,
                                                                       3344, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10439, 0, 3, 9539,
                                                                       2939, 9629, 462, 483,
                                                                       3515, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10565, 0, 3, 9629,
                                                                       2984, 9719, 483, 504,
                                                                       3578, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10691, 0, 3, 9719,
                                                                       3029, 9809, 504, 525,
                                                                       3641, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10817, 0, 3, 9809,
                                                                       3074, 9899, 525, 546,
                                                                       3704, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10943, 0, 3, 9899,
                                                                       3119, 9989, 546, 567,
                                                                       3767, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11069, 0, 3, 9989,
                                                                       3164, 10079, 567, 588,
                                                                       3830, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11195, 0, 3,
                                                                       10079, 3209, 10169, 588,
                                                                       609, 3893, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11321, 0, 3,
                                                                       10169, 3254, 10259, 609,
                                                                       630, 3956, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11447, 0, 3,
                                                                       10259, 3299, 10349, 630,
                                                                       651, 4019, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11573, 0, 3,
                                                                       10439, 3515, 10565, 693,
                                                                       721, 4250, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11741, 0, 3,
                                                                       10565, 3578, 10691, 721,
                                                                       749, 4334, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11909, 0, 3,
                                                                       10691, 3641, 10817, 749,
                                                                       777, 4418, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 12077, 0, 3,
                                                                       10817, 3704, 10943, 777,
                                                                       805, 4502, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 12245, 0, 3,
                                                                       10943, 3767, 11069, 805,
                                                                       833, 4586, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 12413, 0, 3,
                                                                       11069, 3830, 11195, 833,
                                                                       861, 4670, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 12581, 0, 3,
                                                                       11195, 3893, 11321, 861,
                                                                       889, 4754, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 12749, 0, 3,
                                                                       11321, 3956, 11447, 889,
                                                                       917, 4838, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 12917, 0, 3,
                                                                       11573, 4250, 11741, 973,
                                                                       1009, 5138, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 13133, 0, 3,
                                                                       11741, 4334, 11909, 1009,
                                                                       1045, 5246, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 13349, 0, 3,
                                                                       11909, 4418, 12077, 1045,
                                                                       1081, 5354, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 13565, 0, 3,
                                                                       12077, 4502, 12245, 1081,
                                                                       1117, 5462, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 13781, 0, 3,
                                                                       12245, 4586, 12413, 1117,
                                                                       1153, 5570, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 13997, 0, 3,
                                                                       12413, 4670, 12581, 1153,
                                                                       1189, 5678, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 14213, 0, 3,
                                                                       12581, 4754, 12749, 1189,
                                                                       1225, 5786, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 14429, 0, 3,
                                                                       12917, 5138, 13133, 1297,
                                                                       1342, 6164, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 14699, 0, 3,
                                                                       13133, 5246, 13349, 1342,
                                                                       1387, 6299, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 14969, 0, 3,
                                                                       13349, 5354, 13565, 1387,
                                                                       1432, 6434, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 15239, 0, 3,
                                                                       13565, 5462, 13781, 1432,
                                                                       1477, 6569, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 15509, 0, 3,
                                                                       13781, 5570, 13997, 1477,
                                                                       1522, 6704, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 15779, 0, 3,
                                                                       13997, 5678, 14213, 1522,
                                                                       1567, 6839, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 16049, 0, 3,
                                                                       14429, 6164, 14699, 1657,
                                                                       1712, 7304, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 16379, 0, 3,
                                                                       14699, 6299, 14969, 1712,
                                                                       1767, 7469, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 16709, 0, 3,
                                                                       14969, 6434, 15239, 1767,
                                                                       1822, 7634, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 17039, 0, 3,
                                                                       15239, 6569, 15509, 1822,
                                                                       1877, 7799, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 17369, 0, 3,
                                                                       15509, 6704, 15779, 1877,
                                                                       1932, 7964, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17699, 3, 2042,
                                                                       2045, 8129, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17709, 3, 2045,
                                                                       2048, 8135, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17719, 3, 2048,
                                                                       2051, 8141, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17729, 3, 2051,
                                                                       2054, 8147, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17739, 3, 2054,
                                                                       2057, 8153, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17749, 3, 2057,
                                                                       2060, 8159, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17759, 3, 2060,
                                                                       2063, 8165, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17769, 3, 2063,
                                                                       2066, 8171, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17779, 3, 2066,
                                                                       2069, 8177, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17789, 3, 2069,
                                                                       2072, 8183, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17799, 3, 2072,
                                                                       2075, 8189, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17809, 3, 2075,
                                                                       2078, 8195, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17819, 3, 2078,
                                                                       2081, 8201, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17829, 3, 2081,
                                                                       2084, 8207, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 17839, 0, 3,
                                                                       17699, 8129, 17709, 8213,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 17869, 0, 3,
                                                                       17709, 8135, 17719, 8231,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 17899, 0, 3,
                                                                       17719, 8141, 17729, 8249,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 17929, 0, 3,
                                                                       17729, 8147, 17739, 8267,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 17959, 0, 3,
                                                                       17739, 8153, 17749, 8285,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 17989, 0, 3,
                                                                       17749, 8159, 17759, 8303,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18019, 0, 3,
                                                                       17759, 8165, 17769, 8321,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18049, 0, 3,
                                                                       17769, 8171, 17779, 8339,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18079, 0, 3,
                                                                       17779, 8177, 17789, 8357,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18109, 0, 3,
                                                                       17789, 8183, 17799, 8375,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18139, 0, 3,
                                                                       17799, 8189, 17809, 8393,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18169, 0, 3,
                                                                       17809, 8195, 17819, 8411,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18199, 0, 3,
                                                                       17819, 8201, 17829, 8429,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18229, 0, 3,
                                                                       17839, 8213, 17869, 2207,
                                                                       2225, 8447, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18289, 0, 3,
                                                                       17869, 8231, 17899, 2225,
                                                                       2243, 8483, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18349, 0, 3,
                                                                       17899, 8249, 17929, 2243,
                                                                       2261, 8519, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18409, 0, 3,
                                                                       17929, 8267, 17959, 2261,
                                                                       2279, 8555, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18469, 0, 3,
                                                                       17959, 8285, 17989, 2279,
                                                                       2297, 8591, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18529, 0, 3,
                                                                       17989, 8303, 18019, 2297,
                                                                       2315, 8627, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18589, 0, 3,
                                                                       18019, 8321, 18049, 2315,
                                                                       2333, 8663, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18649, 0, 3,
                                                                       18049, 8339, 18079, 2333,
                                                                       2351, 8699, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18709, 0, 3,
                                                                       18079, 8357, 18109, 2351,
                                                                       2369, 8735, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18769, 0, 3,
                                                                       18109, 8375, 18139, 2369,
                                                                       2387, 8771, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18829, 0, 3,
                                                                       18139, 8393, 18169, 2387,
                                                                       2405, 8807, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18889, 0, 3,
                                                                       18169, 8411, 18199, 2405,
                                                                       2423, 8843, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 18949, 0, 3,
                                                                       18229, 8447, 18289, 2459,
                                                                       2489, 8879, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19049, 0, 3,
                                                                       18289, 8483, 18349, 2489,
                                                                       2519, 8939, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19149, 0, 3,
                                                                       18349, 8519, 18409, 2519,
                                                                       2549, 8999, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19249, 0, 3,
                                                                       18409, 8555, 18469, 2549,
                                                                       2579, 9059, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19349, 0, 3,
                                                                       18469, 8591, 18529, 2579,
                                                                       2609, 9119, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19449, 0, 3,
                                                                       18529, 8627, 18589, 2609,
                                                                       2639, 9179, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19549, 0, 3,
                                                                       18589, 8663, 18649, 2639,
                                                                       2669, 9239, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19649, 0, 3,
                                                                       18649, 8699, 18709, 2669,
                                                                       2699, 9299, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19749, 0, 3,
                                                                       18709, 8735, 18769, 2699,
                                                                       2729, 9359, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19849, 0, 3,
                                                                       18769, 8771, 18829, 2729,
                                                                       2759, 9419, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19949, 0, 3,
                                                                       18829, 8807, 18889, 2759,
                                                                       2789, 9479, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 20049, 0, 3,
                                                                       18949, 8879, 19049, 2849,
                                                                       2894, 9539, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 20199, 0, 3,
                                                                       19049, 8939, 19149, 2894,
                                                                       2939, 9629, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 20349, 0, 3,
                                                                       19149, 8999, 19249, 2939,
                                                                       2984, 9719, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 20499, 0, 3,
                                                                       19249, 9059, 19349, 2984,
                                                                       3029, 9809, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 20649, 0, 3,
                                                                       19349, 9119, 19449, 3029,
                                                                       3074, 9899, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 20799, 0, 3,
                                                                       19449, 9179, 19549, 3074,
                                                                       3119, 9989, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 20949, 0, 3,
                                                                       19549, 9239, 19649, 3119,
                                                                       3164, 10079, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 21099, 0, 3,
                                                                       19649, 9299, 19749, 3164,
                                                                       3209, 10169, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 21249, 0, 3,
                                                                       19749, 9359, 19849, 3209,
                                                                       3254, 10259, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 21399, 0, 3,
                                                                       19849, 9419, 19949, 3254,
                                                                       3299, 10349, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 21549, 0, 3,
                                                                       20049, 9539, 20199, 3389,
                                                                       3452, 10439, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 21759, 0, 3,
                                                                       20199, 9629, 20349, 3452,
                                                                       3515, 10565, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 21969, 0, 3,
                                                                       20349, 9719, 20499, 3515,
                                                                       3578, 10691, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 22179, 0, 3,
                                                                       20499, 9809, 20649, 3578,
                                                                       3641, 10817, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 22389, 0, 3,
                                                                       20649, 9899, 20799, 3641,
                                                                       3704, 10943, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 22599, 0, 3,
                                                                       20799, 9989, 20949, 3704,
                                                                       3767, 11069, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 22809, 0, 3,
                                                                       20949, 10079, 21099, 3767,
                                                                       3830, 11195, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 23019, 0, 3,
                                                                       21099, 10169, 21249, 3830,
                                                                       3893, 11321, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 23229, 0, 3,
                                                                       21249, 10259, 21399, 3893,
                                                                       3956, 11447, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 23439, 0, 3,
                                                                       21549, 10439, 21759, 4082,
                                                                       4166, 11573, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 23719, 0, 3,
                                                                       21759, 10565, 21969, 4166,
                                                                       4250, 11741, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 23999, 0, 3,
                                                                       21969, 10691, 22179, 4250,
                                                                       4334, 11909, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 24279, 0, 3,
                                                                       22179, 10817, 22389, 4334,
                                                                       4418, 12077, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 24559, 0, 3,
                                                                       22389, 10943, 22599, 4418,
                                                                       4502, 12245, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 24839, 0, 3,
                                                                       22599, 11069, 22809, 4502,
                                                                       4586, 12413, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 25119, 0, 3,
                                                                       22809, 11195, 23019, 4586,
                                                                       4670, 12581, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 25399, 0, 3,
                                                                       23019, 11321, 23229, 4670,
                                                                       4754, 12749, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 25679, 0, 3,
                                                                       23439, 11573, 23719, 4922,
                                                                       5030, 12917, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 26039, 0, 3,
                                                                       23719, 11741, 23999, 5030,
                                                                       5138, 13133, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 26399, 0, 3,
                                                                       23999, 11909, 24279, 5138,
                                                                       5246, 13349, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 26759, 0, 3,
                                                                       24279, 12077, 24559, 5246,
                                                                       5354, 13565, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 27119, 0, 3,
                                                                       24559, 12245, 24839, 5354,
                                                                       5462, 13781, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 27479, 0, 3,
                                                                       24839, 12413, 25119, 5462,
                                                                       5570, 13997, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 27839, 0, 3,
                                                                       25119, 12581, 25399, 5570,
                                                                       5678, 14213, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 28199, 0, 3,
                                                                       25679, 12917, 26039, 5894,
                                                                       6029, 14429, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 28649, 0, 3,
                                                                       26039, 13133, 26399, 6029,
                                                                       6164, 14699, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 29099, 0, 3,
                                                                       26399, 13349, 26759, 6164,
                                                                       6299, 14969, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 29549, 0, 3,
                                                                       26759, 13565, 27119, 6299,
                                                                       6434, 15239, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 29999, 0, 3,
                                                                       27119, 13781, 27479, 6434,
                                                                       6569, 15509, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 30449, 0, 3,
                                                                       27479, 13997, 27839, 6569,
                                                                       6704, 15779, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 30899, 0, 3,
                                                                       28199, 14429, 28649, 6974,
                                                                       7139, 16049, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 31449, 0, 3,
                                                                       28649, 14699, 29099, 7139,
                                                                       7304, 16379, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 31999, 0, 3,
                                                                       29099, 14969, 29549, 7304,
                                                                       7469, 16709, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 32549, 0, 3,
                                                                       29549, 15239, 29999, 7469,
                                                                       7634, 17039, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 33099, 0, 3,
                                                                       29999, 15509, 30449, 7634,
                                                                       7799, 17369, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33649, 3, 8129,
                                                                       8135, 17719, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33664, 3, 8135,
                                                                       8141, 17729, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33679, 3, 8141,
                                                                       8147, 17739, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33694, 3, 8147,
                                                                       8153, 17749, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33709, 3, 8153,
                                                                       8159, 17759, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33724, 3, 8159,
                                                                       8165, 17769, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33739, 3, 8165,
                                                                       8171, 17779, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33754, 3, 8171,
                                                                       8177, 17789, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33769, 3, 8177,
                                                                       8183, 17799, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33784, 3, 8183,
                                                                       8189, 17809, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33799, 3, 8189,
                                                                       8195, 17819, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33814, 3, 8195,
                                                                       8201, 17829, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 33829, 0, 3,
                                                                       33649, 17719, 33664, 8213,
                                                                       8231, 17899, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 33874, 0, 3,
                                                                       33664, 17729, 33679, 8231,
                                                                       8249, 17929, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 33919, 0, 3,
                                                                       33679, 17739, 33694, 8249,
                                                                       8267, 17959, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 33964, 0, 3,
                                                                       33694, 17749, 33709, 8267,
                                                                       8285, 17989, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 34009, 0, 3,
                                                                       33709, 17759, 33724, 8285,
                                                                       8303, 18019, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 34054, 0, 3,
                                                                       33724, 17769, 33739, 8303,
                                                                       8321, 18049, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 34099, 0, 3,
                                                                       33739, 17779, 33754, 8321,
                                                                       8339, 18079, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 34144, 0, 3,
                                                                       33754, 17789, 33769, 8339,
                                                                       8357, 18109, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 34189, 0, 3,
                                                                       33769, 17799, 33784, 8357,
                                                                       8375, 18139, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 34234, 0, 3,
                                                                       33784, 17809, 33799, 8375,
                                                                       8393, 18169, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 34279, 0, 3,
                                                                       33799, 17819, 33814, 8393,
                                                                       8411, 18199, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 34324, 0, 3,
                                                                       33829, 17899, 33874, 8447,
                                                                       8483, 18349, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 34414, 0, 3,
                                                                       33874, 17929, 33919, 8483,
                                                                       8519, 18409, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 34504, 0, 3,
                                                                       33919, 17959, 33964, 8519,
                                                                       8555, 18469, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 34594, 0, 3,
                                                                       33964, 17989, 34009, 8555,
                                                                       8591, 18529, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 34684, 0, 3,
                                                                       34009, 18019, 34054, 8591,
                                                                       8627, 18589, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 34774, 0, 3,
                                                                       34054, 18049, 34099, 8627,
                                                                       8663, 18649, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 34864, 0, 3,
                                                                       34099, 18079, 34144, 8663,
                                                                       8699, 18709, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 34954, 0, 3,
                                                                       34144, 18109, 34189, 8699,
                                                                       8735, 18769, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 35044, 0, 3,
                                                                       34189, 18139, 34234, 8735,
                                                                       8771, 18829, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 35134, 0, 3,
                                                                       34234, 18169, 34279, 8771,
                                                                       8807, 18889, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 35224, 0, 3,
                                                                       34324, 18349, 34414, 8879,
                                                                       8939, 19149, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 35374, 0, 3,
                                                                       34414, 18409, 34504, 8939,
                                                                       8999, 19249, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 35524, 0, 3,
                                                                       34504, 18469, 34594, 8999,
                                                                       9059, 19349, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 35674, 0, 3,
                                                                       34594, 18529, 34684, 9059,
                                                                       9119, 19449, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 35824, 0, 3,
                                                                       34684, 18589, 34774, 9119,
                                                                       9179, 19549, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 35974, 0, 3,
                                                                       34774, 18649, 34864, 9179,
                                                                       9239, 19649, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 36124, 0, 3,
                                                                       34864, 18709, 34954, 9239,
                                                                       9299, 19749, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 36274, 0, 3,
                                                                       34954, 18769, 35044, 9299,
                                                                       9359, 19849, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 36424, 0, 3,
                                                                       35044, 18829, 35134, 9359,
                                                                       9419, 19949, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 36574, 0, 3,
                                                                       35224, 19149, 35374, 9539,
                                                                       9629, 20349, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 36799, 0, 3,
                                                                       35374, 19249, 35524, 9629,
                                                                       9719, 20499, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 37024, 0, 3,
                                                                       35524, 19349, 35674, 9719,
                                                                       9809, 20649, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 37249, 0, 3,
                                                                       35674, 19449, 35824, 9809,
                                                                       9899, 20799, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 37474, 0, 3,
                                                                       35824, 19549, 35974, 9899,
                                                                       9989, 20949, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 37699, 0, 3,
                                                                       35974, 19649, 36124, 9989,
                                                                       10079, 21099, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 37924, 0, 3,
                                                                       36124, 19749, 36274,
                                                                       10079, 10169, 21249,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 38149, 0, 3,
                                                                       36274, 19849, 36424,
                                                                       10169, 10259, 21399,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 38374, 0, 3,
                                                                       36574, 20349, 36799,
                                                                       10439, 10565, 21969,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 38689, 0, 3,
                                                                       36799, 20499, 37024,
                                                                       10565, 10691, 22179,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 39004, 0, 3,
                                                                       37024, 20649, 37249,
                                                                       10691, 10817, 22389,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 39319, 0, 3,
                                                                       37249, 20799, 37474,
                                                                       10817, 10943, 22599,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 39634, 0, 3,
                                                                       37474, 20949, 37699,
                                                                       10943, 11069, 22809,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 39949, 0, 3,
                                                                       37699, 21099, 37924,
                                                                       11069, 11195, 23019,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 40264, 0, 3,
                                                                       37924, 21249, 38149,
                                                                       11195, 11321, 23229,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 40579, 0, 3,
                                                                       38374, 21969, 38689,
                                                                       11573, 11741, 23999,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 40999, 0, 3,
                                                                       38689, 22179, 39004,
                                                                       11741, 11909, 24279,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 41419, 0, 3,
                                                                       39004, 22389, 39319,
                                                                       11909, 12077, 24559,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 41839, 0, 3,
                                                                       39319, 22599, 39634,
                                                                       12077, 12245, 24839,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 42259, 0, 3,
                                                                       39634, 22809, 39949,
                                                                       12245, 12413, 25119,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 42679, 0, 3,
                                                                       39949, 23019, 40264,
                                                                       12413, 12581, 25399,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 43099, 0, 3,
                                                                       40579, 23999, 40999,
                                                                       12917, 13133, 26399,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 43639, 0, 3,
                                                                       40999, 24279, 41419,
                                                                       13133, 13349, 26759,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 44179, 0, 3,
                                                                       41419, 24559, 41839,
                                                                       13349, 13565, 27119,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 44719, 0, 3,
                                                                       41839, 24839, 42259,
                                                                       13565, 13781, 27479,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 45259, 0, 3,
                                                                       42259, 25119, 42679,
                                                                       13781, 13997, 27839,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 45799, 0, 3,
                                                                       43099, 26399, 43639,
                                                                       14429, 14699, 29099,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 46474, 0, 3,
                                                                       43639, 26759, 44179,
                                                                       14699, 14969, 29549,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 47149, 0, 3,
                                                                       44179, 27119, 44719,
                                                                       14969, 15239, 29999,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 47824, 0, 3,
                                                                       44719, 27479, 45259,
                                                                       15239, 15509, 30449,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 48499, 0, 3,
                                                                       45799, 29099, 46474,
                                                                       16049, 16379, 31999,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 49324, 0, 3,
                                                                       46474, 29549, 47149,
                                                                       16379, 16709, 32549,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 50149, 0, 3,
                                                                       47149, 29999, 47824,
                                                                       16709, 17039, 33099,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 50974, 3, 17699,
                                                                       17709, 33649, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 50995, 3, 17709,
                                                                       17719, 33664, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51016, 3, 17719,
                                                                       17729, 33679, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51037, 3, 17729,
                                                                       17739, 33694, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51058, 3, 17739,
                                                                       17749, 33709, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51079, 3, 17749,
                                                                       17759, 33724, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51100, 3, 17759,
                                                                       17769, 33739, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51121, 3, 17769,
                                                                       17779, 33754, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51142, 3, 17779,
                                                                       17789, 33769, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51163, 3, 17789,
                                                                       17799, 33784, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51184, 3, 17799,
                                                                       17809, 33799, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51205, 3, 17809,
                                                                       17819, 33814, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 51226, 0, 3,
                                                                       50974, 33649, 50995,
                                                                       17839, 17869, 33829,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 51289, 0, 3,
                                                                       50995, 33664, 51016,
                                                                       17869, 17899, 33874,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 51352, 0, 3,
                                                                       51016, 33679, 51037,
                                                                       17899, 17929, 33919,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 51415, 0, 3,
                                                                       51037, 33694, 51058,
                                                                       17929, 17959, 33964,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 51478, 0, 3,
                                                                       51058, 33709, 51079,
                                                                       17959, 17989, 34009,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 51541, 0, 3,
                                                                       51079, 33724, 51100,
                                                                       17989, 18019, 34054,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 51604, 0, 3,
                                                                       51100, 33739, 51121,
                                                                       18019, 18049, 34099,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 51667, 0, 3,
                                                                       51121, 33754, 51142,
                                                                       18049, 18079, 34144,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 51730, 0, 3,
                                                                       51142, 33769, 51163,
                                                                       18079, 18109, 34189,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 51793, 0, 3,
                                                                       51163, 33784, 51184,
                                                                       18109, 18139, 34234,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 51856, 0, 3,
                                                                       51184, 33799, 51205,
                                                                       18139, 18169, 34279,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 51919, 0, 3,
                                                                       51226, 33829, 51289,
                                                                       18229, 18289, 34324,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 52045, 0, 3,
                                                                       51289, 33874, 51352,
                                                                       18289, 18349, 34414,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 52171, 0, 3,
                                                                       51352, 33919, 51415,
                                                                       18349, 18409, 34504,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 52297, 0, 3,
                                                                       51415, 33964, 51478,
                                                                       18409, 18469, 34594,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 52423, 0, 3,
                                                                       51478, 34009, 51541,
                                                                       18469, 18529, 34684,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 52549, 0, 3,
                                                                       51541, 34054, 51604,
                                                                       18529, 18589, 34774,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 52675, 0, 3,
                                                                       51604, 34099, 51667,
                                                                       18589, 18649, 34864,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 52801, 0, 3,
                                                                       51667, 34144, 51730,
                                                                       18649, 18709, 34954,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 52927, 0, 3,
                                                                       51730, 34189, 51793,
                                                                       18709, 18769, 35044,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 53053, 0, 3,
                                                                       51793, 34234, 51856,
                                                                       18769, 18829, 35134,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 53179, 0, 3,
                                                                       51919, 34324, 52045,
                                                                       18949, 19049, 35224,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 53389, 0, 3,
                                                                       52045, 34414, 52171,
                                                                       19049, 19149, 35374,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 53599, 0, 3,
                                                                       52171, 34504, 52297,
                                                                       19149, 19249, 35524,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 53809, 0, 3,
                                                                       52297, 34594, 52423,
                                                                       19249, 19349, 35674,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 54019, 0, 3,
                                                                       52423, 34684, 52549,
                                                                       19349, 19449, 35824,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 54229, 0, 3,
                                                                       52549, 34774, 52675,
                                                                       19449, 19549, 35974,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 54439, 0, 3,
                                                                       52675, 34864, 52801,
                                                                       19549, 19649, 36124,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 54649, 0, 3,
                                                                       52801, 34954, 52927,
                                                                       19649, 19749, 36274,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 54859, 0, 3,
                                                                       52927, 35044, 53053,
                                                                       19749, 19849, 36424,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 55069, 0, 3,
                                                                       53179, 35224, 53389,
                                                                       20049, 20199, 36574,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 55384, 0, 3,
                                                                       53389, 35374, 53599,
                                                                       20199, 20349, 36799,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 55699, 0, 3,
                                                                       53599, 35524, 53809,
                                                                       20349, 20499, 37024,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 56014, 0, 3,
                                                                       53809, 35674, 54019,
                                                                       20499, 20649, 37249,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 56329, 0, 3,
                                                                       54019, 35824, 54229,
                                                                       20649, 20799, 37474,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 56644, 0, 3,
                                                                       54229, 35974, 54439,
                                                                       20799, 20949, 37699,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 56959, 0, 3,
                                                                       54439, 36124, 54649,
                                                                       20949, 21099, 37924,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 57274, 0, 3,
                                                                       54649, 36274, 54859,
                                                                       21099, 21249, 38149,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 57589, 0, 3,
                                                                       55069, 36574, 55384,
                                                                       21549, 21759, 38374,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 58030, 0, 3,
                                                                       55384, 36799, 55699,
                                                                       21759, 21969, 38689,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 58471, 0, 3,
                                                                       55699, 37024, 56014,
                                                                       21969, 22179, 39004,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 58912, 0, 3,
                                                                       56014, 37249, 56329,
                                                                       22179, 22389, 39319,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 59353, 0, 3,
                                                                       56329, 37474, 56644,
                                                                       22389, 22599, 39634,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 59794, 0, 3,
                                                                       56644, 37699, 56959,
                                                                       22599, 22809, 39949,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 60235, 0, 3,
                                                                       56959, 37924, 57274,
                                                                       22809, 23019, 40264,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 60676, 0, 3,
                                                                       57589, 38374, 58030,
                                                                       23439, 23719, 40579,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 61264, 0, 3,
                                                                       58030, 38689, 58471,
                                                                       23719, 23999, 40999,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 61852, 0, 3,
                                                                       58471, 39004, 58912,
                                                                       23999, 24279, 41419,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 62440, 0, 3,
                                                                       58912, 39319, 59353,
                                                                       24279, 24559, 41839,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 63028, 0, 3,
                                                                       59353, 39634, 59794,
                                                                       24559, 24839, 42259,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 63616, 0, 3,
                                                                       59794, 39949, 60235,
                                                                       24839, 25119, 42679,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 64204, 0, 3,
                                                                       60676, 40579, 61264,
                                                                       25679, 26039, 43099,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 64960, 0, 3,
                                                                       61264, 40999, 61852,
                                                                       26039, 26399, 43639,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 65716, 0, 3,
                                                                       61852, 41419, 62440,
                                                                       26399, 26759, 44179,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 66472, 0, 3,
                                                                       62440, 41839, 63028,
                                                                       26759, 27119, 44719,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 67228, 0, 3,
                                                                       63028, 42259, 63616,
                                                                       27119, 27479, 45259,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 67984, 0, 3,
                                                                       64204, 43099, 64960,
                                                                       28199, 28649, 45799,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 68929, 0, 3,
                                                                       64960, 43639, 65716,
                                                                       28649, 29099, 46474,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 69874, 0, 3,
                                                                       65716, 44179, 66472,
                                                                       29099, 29549, 47149,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 70819, 0, 3,
                                                                       66472, 44719, 67228,
                                                                       29549, 29999, 47824,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 71764, 0, 3,
                                                                       67984, 45799, 68929,
                                                                       30899, 31449, 48499,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 72919, 0, 3,
                                                                       68929, 46474, 69874,
                                                                       31449, 31999, 49324,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 74074, 0, 3,
                                                                       69874, 47149, 70819,
                                                                       31999, 32549, 50149,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75229, 3, 33649,
                                                                       33664, 51016, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75257, 3, 33664,
                                                                       33679, 51037, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75285, 3, 33679,
                                                                       33694, 51058, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75313, 3, 33694,
                                                                       33709, 51079, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75341, 3, 33709,
                                                                       33724, 51100, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75369, 3, 33724,
                                                                       33739, 51121, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75397, 3, 33739,
                                                                       33754, 51142, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75425, 3, 33754,
                                                                       33769, 51163, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75453, 3, 33769,
                                                                       33784, 51184, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75481, 3, 33784,
                                                                       33799, 51205, ncols,
                                                                       gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 75509, 0, 3,
                                                                       75229, 51016, 75257,
                                                                       33829, 33874, 51352,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 75593, 0, 3,
                                                                       75257, 51037, 75285,
                                                                       33874, 33919, 51415,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 75677, 0, 3,
                                                                       75285, 51058, 75313,
                                                                       33919, 33964, 51478,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 75761, 0, 3,
                                                                       75313, 51079, 75341,
                                                                       33964, 34009, 51541,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 75845, 0, 3,
                                                                       75341, 51100, 75369,
                                                                       34009, 34054, 51604,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 75929, 0, 3,
                                                                       75369, 51121, 75397,
                                                                       34054, 34099, 51667,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 76013, 0, 3,
                                                                       75397, 51142, 75425,
                                                                       34099, 34144, 51730,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 76097, 0, 3,
                                                                       75425, 51163, 75453,
                                                                       34144, 34189, 51793,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 76181, 0, 3,
                                                                       75453, 51184, 75481,
                                                                       34189, 34234, 51856,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 76265, 0, 3,
                                                                       75509, 51352, 75593,
                                                                       34324, 34414, 52171,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 76433, 0, 3,
                                                                       75593, 51415, 75677,
                                                                       34414, 34504, 52297,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 76601, 0, 3,
                                                                       75677, 51478, 75761,
                                                                       34504, 34594, 52423,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 76769, 0, 3,
                                                                       75761, 51541, 75845,
                                                                       34594, 34684, 52549,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 76937, 0, 3,
                                                                       75845, 51604, 75929,
                                                                       34684, 34774, 52675,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 77105, 0, 3,
                                                                       75929, 51667, 76013,
                                                                       34774, 34864, 52801,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 77273, 0, 3,
                                                                       76013, 51730, 76097,
                                                                       34864, 34954, 52927,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 77441, 0, 3,
                                                                       76097, 51793, 76181,
                                                                       34954, 35044, 53053,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 77609, 0, 3,
                                                                       76265, 52171, 76433,
                                                                       35224, 35374, 53599,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 77889, 0, 3,
                                                                       76433, 52297, 76601,
                                                                       35374, 35524, 53809,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 78169, 0, 3,
                                                                       76601, 52423, 76769,
                                                                       35524, 35674, 54019,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 78449, 0, 3,
                                                                       76769, 52549, 76937,
                                                                       35674, 35824, 54229,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 78729, 0, 3,
                                                                       76937, 52675, 77105,
                                                                       35824, 35974, 54439,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 79009, 0, 3,
                                                                       77105, 52801, 77273,
                                                                       35974, 36124, 54649,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 79289, 0, 3,
                                                                       77273, 52927, 77441,
                                                                       36124, 36274, 54859,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 79569, 0, 3,
                                                                       77609, 53599, 77889,
                                                                       36574, 36799, 55699,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 79989, 0, 3,
                                                                       77889, 53809, 78169,
                                                                       36799, 37024, 56014,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 80409, 0, 3,
                                                                       78169, 54019, 78449,
                                                                       37024, 37249, 56329,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 80829, 0, 3,
                                                                       78449, 54229, 78729,
                                                                       37249, 37474, 56644,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 81249, 0, 3,
                                                                       78729, 54439, 79009,
                                                                       37474, 37699, 56959,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 81669, 0, 3,
                                                                       79009, 54649, 79289,
                                                                       37699, 37924, 57274,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 82089, 0, 3,
                                                                       79569, 55699, 79989,
                                                                       38374, 38689, 58471,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 82677, 0, 3,
                                                                       79989, 56014, 80409,
                                                                       38689, 39004, 58912,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 83265, 0, 3,
                                                                       80409, 56329, 80829,
                                                                       39004, 39319, 59353,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 83853, 0, 3,
                                                                       80829, 56644, 81249,
                                                                       39319, 39634, 59794,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 84441, 0, 3,
                                                                       81249, 56959, 81669,
                                                                       39634, 39949, 60235,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 85029, 0, 3,
                                                                       82089, 58471, 82677,
                                                                       40579, 40999, 61852,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 85813, 0, 3,
                                                                       82677, 58912, 83265,
                                                                       40999, 41419, 62440,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 86597, 0, 3,
                                                                       83265, 59353, 83853,
                                                                       41419, 41839, 63028,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 87381, 0, 3,
                                                                       83853, 59794, 84441,
                                                                       41839, 42259, 63616,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 88165, 0, 3,
                                                                       85029, 61852, 85813,
                                                                       43099, 43639, 65716,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 89173, 0, 3,
                                                                       85813, 62440, 86597,
                                                                       43639, 44179, 66472,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 90181, 0, 3,
                                                                       86597, 63028, 87381,
                                                                       44179, 44719, 67228,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 91189, 0, 3,
                                                                       88165, 65716, 89173,
                                                                       45799, 46474, 69874,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 92449, 0, 3,
                                                                       89173, 66472, 90181,
                                                                       46474, 47149, 70819,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 93709, 0, 3,
                                                                       91189, 69874, 92449,
                                                                       48499, 49324, 74074,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 95249, 3, 50974,
                                                                       50995, 75229, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 95285, 3, 50995,
                                                                       51016, 75257, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 95321, 3, 51016,
                                                                       51037, 75285, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 95357, 3, 51037,
                                                                       51058, 75313, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 95393, 3, 51058,
                                                                       51079, 75341, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 95429, 3, 51079,
                                                                       51100, 75369, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 95465, 3, 51100,
                                                                       51121, 75397, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 95501, 3, 51121,
                                                                       51142, 75425, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 95537, 3, 51142,
                                                                       51163, 75453, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 95573, 3, 51163,
                                                                       51184, 75481, ncols,
                                                                       gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 95609, 0, 3,
                                                                       95249, 75229, 95285,
                                                                       51226, 51289, 75509,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 95717, 0, 3,
                                                                       95285, 75257, 95321,
                                                                       51289, 51352, 75593,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 95825, 0, 3,
                                                                       95321, 75285, 95357,
                                                                       51352, 51415, 75677,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 95933, 0, 3,
                                                                       95357, 75313, 95393,
                                                                       51415, 51478, 75761,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 96041, 0, 3,
                                                                       95393, 75341, 95429,
                                                                       51478, 51541, 75845,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 96149, 0, 3,
                                                                       95429, 75369, 95465,
                                                                       51541, 51604, 75929,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 96257, 0, 3,
                                                                       95465, 75397, 95501,
                                                                       51604, 51667, 76013,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 96365, 0, 3,
                                                                       95501, 75425, 95537,
                                                                       51667, 51730, 76097,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 96473, 0, 3,
                                                                       95537, 75453, 95573,
                                                                       51730, 51793, 76181,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 96581, 0, 3,
                                                                       95609, 75509, 95717,
                                                                       51919, 52045, 76265,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 96797, 0, 3,
                                                                       95717, 75593, 95825,
                                                                       52045, 52171, 76433,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 97013, 0, 3,
                                                                       95825, 75677, 95933,
                                                                       52171, 52297, 76601,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 97229, 0, 3,
                                                                       95933, 75761, 96041,
                                                                       52297, 52423, 76769,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 97445, 0, 3,
                                                                       96041, 75845, 96149,
                                                                       52423, 52549, 76937,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 97661, 0, 3,
                                                                       96149, 75929, 96257,
                                                                       52549, 52675, 77105,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 97877, 0, 3,
                                                                       96257, 76013, 96365,
                                                                       52675, 52801, 77273,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 98093, 0, 3,
                                                                       96365, 76097, 96473,
                                                                       52801, 52927, 77441,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 98309, 0, 3,
                                                                       96581, 76265, 96797,
                                                                       53179, 53389, 77609,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 98669, 0, 3,
                                                                       96797, 76433, 97013,
                                                                       53389, 53599, 77889,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 99029, 0, 3,
                                                                       97013, 76601, 97229,
                                                                       53599, 53809, 78169,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 99389, 0, 3,
                                                                       97229, 76769, 97445,
                                                                       53809, 54019, 78449,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 99749, 0, 3,
                                                                       97445, 76937, 97661,
                                                                       54019, 54229, 78729,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 100109, 0, 3,
                                                                       97661, 77105, 97877,
                                                                       54229, 54439, 79009,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 100469, 0, 3,
                                                                       97877, 77273, 98093,
                                                                       54439, 54649, 79289,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 100829, 0, 3,
                                                                       98309, 77609, 98669,
                                                                       55069, 55384, 79569,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 101369, 0, 3,
                                                                       98669, 77889, 99029,
                                                                       55384, 55699, 79989,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 101909, 0, 3,
                                                                       99029, 78169, 99389,
                                                                       55699, 56014, 80409,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 102449, 0, 3,
                                                                       99389, 78449, 99749,
                                                                       56014, 56329, 80829,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 102989, 0, 3,
                                                                       99749, 78729, 100109,
                                                                       56329, 56644, 81249,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 103529, 0, 3,
                                                                       100109, 79009, 100469,
                                                                       56644, 56959, 81669,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 104069, 0, 3,
                                                                       100829, 79569, 101369,
                                                                       57589, 58030, 82089,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 104825, 0, 3,
                                                                       101369, 79989, 101909,
                                                                       58030, 58471, 82677,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 105581, 0, 3,
                                                                       101909, 80409, 102449,
                                                                       58471, 58912, 83265,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 106337, 0, 3,
                                                                       102449, 80829, 102989,
                                                                       58912, 59353, 83853,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 107093, 0, 3,
                                                                       102989, 81249, 103529,
                                                                       59353, 59794, 84441,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 107849, 0, 3,
                                                                       104069, 82089, 104825,
                                                                       60676, 61264, 85029,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 108857, 0, 3,
                                                                       104825, 82677, 105581,
                                                                       61264, 61852, 85813,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 109865, 0, 3,
                                                                       105581, 83265, 106337,
                                                                       61852, 62440, 86597,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 110873, 0, 3,
                                                                       106337, 83853, 107093,
                                                                       62440, 63028, 87381,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 111881, 0, 3,
                                                                       107849, 85029, 108857,
                                                                       64204, 64960, 88165,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 113177, 0, 3,
                                                                       108857, 85813, 109865,
                                                                       64960, 65716, 89173,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 114473, 0, 3,
                                                                       109865, 86597, 110873,
                                                                       65716, 66472, 90181,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 115769, 0, 3,
                                                                       111881, 88165, 113177,
                                                                       67984, 68929, 91189,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 117389, 0, 3,
                                                                       113177, 89173, 114473,
                                                                       68929, 69874, 92449,
                                                                       ncols, gamma, p, q);

                    compute_prim_smk_three_center_electron_repulsion_0(buffer, 119009, 0, 3,
                                                                       115769, 91189, 117389,
                                                                       71764, 72919, 93709,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 120989, 104069, 756, ncols);

                    simdfunc::contract_primitives(buffer, 122060, 107849, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 123488, 111881, 1296, ncols);

                    simdfunc::contract_primitives(buffer, 125324, 115769, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 127619, 119009, 1980, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 121745, 120989, 21, 1, nmax);

        simdtrf::transform_k_inner(buffer, 123068, 122060, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 124784, 123488, 36, 1, nmax);

        simdtrf::transform_k_inner(buffer, 126944, 125324, 45, 1, nmax);

        simdtrf::transform_k_inner(buffer, 129599, 127619, 55, 1, nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 130424, 121745, 123068, 15, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 131369, 123068, 124784, 15, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 132629, 124784, 126944, 15, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 134249, 126944, 129599, 15, nmax);

        simdtrf::compute_hrr_dh(buffer, coordinates, 136274, 130424, 131369, 15, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 138164, 131369, 132629, 15, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 140684, 132629, 134249, 15, nmax);

        simdtrf::compute_hrr_fh(buffer, coordinates, 143924, 136274, 138164, 15, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 147074, 138164, 140684, 15, nmax);

        simdtrf::compute_hrr_gh(buffer, coordinates, 151274, 143924, 147074, 15, nmax);

        simdtrf::transform_h_inner(buffer, 155999, 151274, 15, 15, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 155999, 165, nmax);
    }

    for (size_t m = 0; m < 1485; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
