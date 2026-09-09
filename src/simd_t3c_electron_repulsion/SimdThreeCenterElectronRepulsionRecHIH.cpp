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


#include "SimdThreeCenterElectronRepulsionRecHIH.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSOD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOS.hpp"
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
#include "SimdTransferDM.hpp"
#include "SimdTransferFI.hpp"
#include "SimdTransferFK.hpp"
#include "SimdTransferFL.hpp"
#include "SimdTransferGI.hpp"
#include "SimdTransferGK.hpp"
#include "SimdTransferHI.hpp"
#include "SimdTransferPI.hpp"
#include "SimdTransferPK.hpp"
#include "SimdTransferPL.hpp"
#include "SimdTransferPM.hpp"
#include "SimdTransferPN.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_hih_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_hih_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 154192, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1573 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 154192, 93901, 8998, dimensions);

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

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2042, 0, 3, 1297,
                                                                       1342, 1657, 1712, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2108, 0, 3, 1342,
                                                                       1387, 1712, 1767, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2174, 0, 3, 1387,
                                                                       1432, 1767, 1822, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2240, 0, 3, 1432,
                                                                       1477, 1822, 1877, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2306, 0, 3, 1477,
                                                                       1522, 1877, 1932, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2372, 0, 3, 1522,
                                                                       1567, 1932, 1987, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 2438, 0, 3, 1657,
                                                                       1712, 2042, 2108, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 2516, 0, 3, 1712,
                                                                       1767, 2108, 2174, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 2594, 0, 3, 1767,
                                                                       1822, 2174, 2240, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 2672, 0, 3, 1822,
                                                                       1877, 2240, 2306, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 2750, 0, 3, 1877,
                                                                       1932, 2306, 2372, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2828, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2831, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2834, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2837, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2840, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2843, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2846, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2849, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2852, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2855, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2858, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2861, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2864, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2867, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2870, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2873, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2876, 3, 9, 29,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2885, 3, 10, 32,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2894, 3, 11, 35,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2903, 3, 12, 38,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2912, 3, 13, 41,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2921, 3, 14, 44,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2930, 3, 15, 47,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2939, 3, 16, 50,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2948, 3, 17, 53,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2957, 3, 18, 56,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2966, 3, 19, 59,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2975, 3, 20, 62,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2984, 3, 21, 65,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2993, 3, 23, 68,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3011, 3, 26, 74,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3029, 3, 29, 80,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3047, 3, 32, 86,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3065, 3, 35, 92,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3083, 3, 38, 98,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3101, 3, 41, 104,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3119, 3, 44, 110,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3137, 3, 47, 116,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3155, 3, 50, 122,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3173, 3, 53, 128,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3191, 3, 56, 134,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3209, 3, 59, 140,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3227, 3, 62, 146,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3245, 3, 68, 152,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3275, 3, 74, 162,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3305, 3, 80, 172,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3335, 3, 86, 182,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3365, 3, 92, 192,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3395, 3, 98, 202,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3425, 3, 104, 212,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3455, 3, 110, 222,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3485, 3, 116, 232,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3515, 3, 122, 242,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3545, 3, 128, 252,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3575, 3, 134, 262,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3605, 3, 140, 272,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3635, 3, 152, 282,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3680, 3, 162, 297,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3725, 3, 172, 312,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3770, 3, 182, 327,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3815, 3, 192, 342,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3860, 3, 202, 357,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3905, 3, 212, 372,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3950, 3, 222, 387,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3995, 3, 232, 402,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4040, 3, 242, 417,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4085, 3, 252, 432,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4130, 3, 262, 447,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4175, 3, 282, 462,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4238, 3, 297, 483,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4301, 3, 312, 504,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4364, 3, 327, 525,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4427, 3, 342, 546,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4490, 3, 357, 567,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4553, 3, 372, 588,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4616, 3, 387, 609,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4679, 3, 402, 630,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4742, 3, 417, 651,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4805, 3, 432, 672,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4868, 3, 462, 693,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4952, 3, 483, 721,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5036, 3, 504, 749,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5120, 3, 525, 777,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5204, 3, 546, 805,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5288, 3, 567, 833,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5372, 3, 588, 861,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5456, 3, 609, 889,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5540, 3, 630, 917,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5624, 3, 651, 945,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5708, 3, 693, 973,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5816, 3, 721,
                                                                       1009, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5924, 3, 749,
                                                                       1045, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6032, 3, 777,
                                                                       1081, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6140, 3, 805,
                                                                       1117, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6248, 3, 833,
                                                                       1153, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6356, 3, 861,
                                                                       1189, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6464, 3, 889,
                                                                       1225, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6572, 3, 917,
                                                                       1261, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6680, 3, 973,
                                                                       1297, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6815, 3, 1009,
                                                                       1342, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6950, 3, 1045,
                                                                       1387, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 7085, 3, 1081,
                                                                       1432, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 7220, 3, 1117,
                                                                       1477, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 7355, 3, 1153,
                                                                       1522, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 7490, 3, 1189,
                                                                       1567, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 7625, 3, 1225,
                                                                       1612, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 7760, 3, 1297,
                                                                       1657, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 7925, 3, 1342,
                                                                       1712, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 8090, 3, 1387,
                                                                       1767, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 8255, 3, 1432,
                                                                       1822, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 8420, 3, 1477,
                                                                       1877, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 8585, 3, 1522,
                                                                       1932, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 8750, 3, 1567,
                                                                       1987, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 8915, 3, 1657,
                                                                       2042, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 9113, 3, 1712,
                                                                       2108, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 9311, 3, 1767,
                                                                       2174, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 9509, 3, 1822,
                                                                       2240, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 9707, 3, 1877,
                                                                       2306, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 9905, 3, 1932,
                                                                       2372, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 10103, 3, 2042,
                                                                       2438, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 10337, 3, 2108,
                                                                       2516, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 10571, 3, 2174,
                                                                       2594, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 10805, 3, 2240,
                                                                       2672, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 11039, 3, 2306,
                                                                       2750, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11273, 3, 7, 8,
                                                                       2834, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11279, 3, 8, 9,
                                                                       2837, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11285, 3, 9, 10,
                                                                       2840, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11291, 3, 10, 11,
                                                                       2843, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11297, 3, 11, 12,
                                                                       2846, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11303, 3, 12, 13,
                                                                       2849, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11309, 3, 13, 14,
                                                                       2852, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11315, 3, 14, 15,
                                                                       2855, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11321, 3, 15, 16,
                                                                       2858, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11327, 3, 16, 17,
                                                                       2861, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11333, 3, 17, 18,
                                                                       2864, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11339, 3, 18, 19,
                                                                       2867, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11345, 3, 19, 20,
                                                                       2870, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11351, 3, 20, 21,
                                                                       2873, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11357, 0, 3,
                                                                       11273, 2834, 11279, 2876,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11375, 0, 3,
                                                                       11279, 2837, 11285, 2885,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11393, 0, 3,
                                                                       11285, 2840, 11291, 2894,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11411, 0, 3,
                                                                       11291, 2843, 11297, 2903,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11429, 0, 3,
                                                                       11297, 2846, 11303, 2912,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11447, 0, 3,
                                                                       11303, 2849, 11309, 2921,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11465, 0, 3,
                                                                       11309, 2852, 11315, 2930,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11483, 0, 3,
                                                                       11315, 2855, 11321, 2939,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11501, 0, 3,
                                                                       11321, 2858, 11327, 2948,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11519, 0, 3,
                                                                       11327, 2861, 11333, 2957,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11537, 0, 3,
                                                                       11333, 2864, 11339, 2966,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11555, 0, 3,
                                                                       11339, 2867, 11345, 2975,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11573, 0, 3,
                                                                       11345, 2870, 11351, 2984,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11591, 0, 3,
                                                                       11357, 2876, 11375, 68,
                                                                       74, 3029, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11627, 0, 3,
                                                                       11375, 2885, 11393, 74,
                                                                       80, 3047, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11663, 0, 3,
                                                                       11393, 2894, 11411, 80,
                                                                       86, 3065, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11699, 0, 3,
                                                                       11411, 2903, 11429, 86,
                                                                       92, 3083, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11735, 0, 3,
                                                                       11429, 2912, 11447, 92,
                                                                       98, 3101, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11771, 0, 3,
                                                                       11447, 2921, 11465, 98,
                                                                       104, 3119, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11807, 0, 3,
                                                                       11465, 2930, 11483, 104,
                                                                       110, 3137, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11843, 0, 3,
                                                                       11483, 2939, 11501, 110,
                                                                       116, 3155, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11879, 0, 3,
                                                                       11501, 2948, 11519, 116,
                                                                       122, 3173, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11915, 0, 3,
                                                                       11519, 2957, 11537, 122,
                                                                       128, 3191, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11951, 0, 3,
                                                                       11537, 2966, 11555, 128,
                                                                       134, 3209, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11987, 0, 3,
                                                                       11555, 2975, 11573, 134,
                                                                       140, 3227, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12023, 0, 3,
                                                                       11591, 3029, 11627, 152,
                                                                       162, 3305, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12083, 0, 3,
                                                                       11627, 3047, 11663, 162,
                                                                       172, 3335, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12143, 0, 3,
                                                                       11663, 3065, 11699, 172,
                                                                       182, 3365, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12203, 0, 3,
                                                                       11699, 3083, 11735, 182,
                                                                       192, 3395, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12263, 0, 3,
                                                                       11735, 3101, 11771, 192,
                                                                       202, 3425, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12323, 0, 3,
                                                                       11771, 3119, 11807, 202,
                                                                       212, 3455, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12383, 0, 3,
                                                                       11807, 3137, 11843, 212,
                                                                       222, 3485, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12443, 0, 3,
                                                                       11843, 3155, 11879, 222,
                                                                       232, 3515, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12503, 0, 3,
                                                                       11879, 3173, 11915, 232,
                                                                       242, 3545, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12563, 0, 3,
                                                                       11915, 3191, 11951, 242,
                                                                       252, 3575, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12623, 0, 3,
                                                                       11951, 3209, 11987, 252,
                                                                       262, 3605, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12683, 0, 3,
                                                                       12023, 3305, 12083, 282,
                                                                       297, 3725, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12773, 0, 3,
                                                                       12083, 3335, 12143, 297,
                                                                       312, 3770, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12863, 0, 3,
                                                                       12143, 3365, 12203, 312,
                                                                       327, 3815, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12953, 0, 3,
                                                                       12203, 3395, 12263, 327,
                                                                       342, 3860, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 13043, 0, 3,
                                                                       12263, 3425, 12323, 342,
                                                                       357, 3905, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 13133, 0, 3,
                                                                       12323, 3455, 12383, 357,
                                                                       372, 3950, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 13223, 0, 3,
                                                                       12383, 3485, 12443, 372,
                                                                       387, 3995, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 13313, 0, 3,
                                                                       12443, 3515, 12503, 387,
                                                                       402, 4040, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 13403, 0, 3,
                                                                       12503, 3545, 12563, 402,
                                                                       417, 4085, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 13493, 0, 3,
                                                                       12563, 3575, 12623, 417,
                                                                       432, 4130, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 13583, 0, 3,
                                                                       12683, 3725, 12773, 462,
                                                                       483, 4301, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 13709, 0, 3,
                                                                       12773, 3770, 12863, 483,
                                                                       504, 4364, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 13835, 0, 3,
                                                                       12863, 3815, 12953, 504,
                                                                       525, 4427, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 13961, 0, 3,
                                                                       12953, 3860, 13043, 525,
                                                                       546, 4490, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 14087, 0, 3,
                                                                       13043, 3905, 13133, 546,
                                                                       567, 4553, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 14213, 0, 3,
                                                                       13133, 3950, 13223, 567,
                                                                       588, 4616, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 14339, 0, 3,
                                                                       13223, 3995, 13313, 588,
                                                                       609, 4679, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 14465, 0, 3,
                                                                       13313, 4040, 13403, 609,
                                                                       630, 4742, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 14591, 0, 3,
                                                                       13403, 4085, 13493, 630,
                                                                       651, 4805, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 14717, 0, 3,
                                                                       13583, 4301, 13709, 693,
                                                                       721, 5036, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 14885, 0, 3,
                                                                       13709, 4364, 13835, 721,
                                                                       749, 5120, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 15053, 0, 3,
                                                                       13835, 4427, 13961, 749,
                                                                       777, 5204, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 15221, 0, 3,
                                                                       13961, 4490, 14087, 777,
                                                                       805, 5288, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 15389, 0, 3,
                                                                       14087, 4553, 14213, 805,
                                                                       833, 5372, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 15557, 0, 3,
                                                                       14213, 4616, 14339, 833,
                                                                       861, 5456, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 15725, 0, 3,
                                                                       14339, 4679, 14465, 861,
                                                                       889, 5540, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 15893, 0, 3,
                                                                       14465, 4742, 14591, 889,
                                                                       917, 5624, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 16061, 0, 3,
                                                                       14717, 5036, 14885, 973,
                                                                       1009, 5924, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 16277, 0, 3,
                                                                       14885, 5120, 15053, 1009,
                                                                       1045, 6032, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 16493, 0, 3,
                                                                       15053, 5204, 15221, 1045,
                                                                       1081, 6140, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 16709, 0, 3,
                                                                       15221, 5288, 15389, 1081,
                                                                       1117, 6248, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 16925, 0, 3,
                                                                       15389, 5372, 15557, 1117,
                                                                       1153, 6356, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 17141, 0, 3,
                                                                       15557, 5456, 15725, 1153,
                                                                       1189, 6464, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 17357, 0, 3,
                                                                       15725, 5540, 15893, 1189,
                                                                       1225, 6572, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 17573, 0, 3,
                                                                       16061, 5924, 16277, 1297,
                                                                       1342, 6950, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 17843, 0, 3,
                                                                       16277, 6032, 16493, 1342,
                                                                       1387, 7085, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 18113, 0, 3,
                                                                       16493, 6140, 16709, 1387,
                                                                       1432, 7220, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 18383, 0, 3,
                                                                       16709, 6248, 16925, 1432,
                                                                       1477, 7355, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 18653, 0, 3,
                                                                       16925, 6356, 17141, 1477,
                                                                       1522, 7490, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 18923, 0, 3,
                                                                       17141, 6464, 17357, 1522,
                                                                       1567, 7625, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 19193, 0, 3,
                                                                       17573, 6950, 17843, 1657,
                                                                       1712, 8090, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 19523, 0, 3,
                                                                       17843, 7085, 18113, 1712,
                                                                       1767, 8255, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 19853, 0, 3,
                                                                       18113, 7220, 18383, 1767,
                                                                       1822, 8420, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 20183, 0, 3,
                                                                       18383, 7355, 18653, 1822,
                                                                       1877, 8585, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 20513, 0, 3,
                                                                       18653, 7490, 18923, 1877,
                                                                       1932, 8750, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 20843, 0, 3,
                                                                       19193, 8090, 19523, 2042,
                                                                       2108, 9311, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 21239, 0, 3,
                                                                       19523, 8255, 19853, 2108,
                                                                       2174, 9509, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 21635, 0, 3,
                                                                       19853, 8420, 20183, 2174,
                                                                       2240, 9707, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 22031, 0, 3,
                                                                       20183, 8585, 20513, 2240,
                                                                       2306, 9905, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 22427, 0, 3,
                                                                       20843, 9311, 21239, 2438,
                                                                       2516, 10571, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 22895, 0, 3,
                                                                       21239, 9509, 21635, 2516,
                                                                       2594, 10805, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 23363, 0, 3,
                                                                       21635, 9707, 22031, 2594,
                                                                       2672, 11039, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23831, 3, 2828,
                                                                       2831, 11273, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23841, 3, 2831,
                                                                       2834, 11279, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23851, 3, 2834,
                                                                       2837, 11285, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23861, 3, 2837,
                                                                       2840, 11291, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23871, 3, 2840,
                                                                       2843, 11297, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23881, 3, 2843,
                                                                       2846, 11303, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23891, 3, 2846,
                                                                       2849, 11309, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23901, 3, 2849,
                                                                       2852, 11315, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23911, 3, 2852,
                                                                       2855, 11321, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23921, 3, 2855,
                                                                       2858, 11327, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23931, 3, 2858,
                                                                       2861, 11333, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23941, 3, 2861,
                                                                       2864, 11339, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23951, 3, 2864,
                                                                       2867, 11345, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23961, 3, 2867,
                                                                       2870, 11351, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 23971, 0, 3,
                                                                       23831, 11273, 23841,
                                                                       11357, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24001, 0, 3,
                                                                       23841, 11279, 23851,
                                                                       11375, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24031, 0, 3,
                                                                       23851, 11285, 23861,
                                                                       11393, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24061, 0, 3,
                                                                       23861, 11291, 23871,
                                                                       11411, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24091, 0, 3,
                                                                       23871, 11297, 23881,
                                                                       11429, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24121, 0, 3,
                                                                       23881, 11303, 23891,
                                                                       11447, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24151, 0, 3,
                                                                       23891, 11309, 23901,
                                                                       11465, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24181, 0, 3,
                                                                       23901, 11315, 23911,
                                                                       11483, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24211, 0, 3,
                                                                       23911, 11321, 23921,
                                                                       11501, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24241, 0, 3,
                                                                       23921, 11327, 23931,
                                                                       11519, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24271, 0, 3,
                                                                       23931, 11333, 23941,
                                                                       11537, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24301, 0, 3,
                                                                       23941, 11339, 23951,
                                                                       11555, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24331, 0, 3,
                                                                       23951, 11345, 23961,
                                                                       11573, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24361, 0, 3,
                                                                       23971, 11357, 24001, 2993,
                                                                       3011, 11591, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24421, 0, 3,
                                                                       24001, 11375, 24031, 3011,
                                                                       3029, 11627, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24481, 0, 3,
                                                                       24031, 11393, 24061, 3029,
                                                                       3047, 11663, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24541, 0, 3,
                                                                       24061, 11411, 24091, 3047,
                                                                       3065, 11699, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24601, 0, 3,
                                                                       24091, 11429, 24121, 3065,
                                                                       3083, 11735, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24661, 0, 3,
                                                                       24121, 11447, 24151, 3083,
                                                                       3101, 11771, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24721, 0, 3,
                                                                       24151, 11465, 24181, 3101,
                                                                       3119, 11807, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24781, 0, 3,
                                                                       24181, 11483, 24211, 3119,
                                                                       3137, 11843, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24841, 0, 3,
                                                                       24211, 11501, 24241, 3137,
                                                                       3155, 11879, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24901, 0, 3,
                                                                       24241, 11519, 24271, 3155,
                                                                       3173, 11915, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24961, 0, 3,
                                                                       24271, 11537, 24301, 3173,
                                                                       3191, 11951, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 25021, 0, 3,
                                                                       24301, 11555, 24331, 3191,
                                                                       3209, 11987, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25081, 0, 3,
                                                                       24361, 11591, 24421, 3245,
                                                                       3275, 12023, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25181, 0, 3,
                                                                       24421, 11627, 24481, 3275,
                                                                       3305, 12083, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25281, 0, 3,
                                                                       24481, 11663, 24541, 3305,
                                                                       3335, 12143, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25381, 0, 3,
                                                                       24541, 11699, 24601, 3335,
                                                                       3365, 12203, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25481, 0, 3,
                                                                       24601, 11735, 24661, 3365,
                                                                       3395, 12263, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25581, 0, 3,
                                                                       24661, 11771, 24721, 3395,
                                                                       3425, 12323, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25681, 0, 3,
                                                                       24721, 11807, 24781, 3425,
                                                                       3455, 12383, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25781, 0, 3,
                                                                       24781, 11843, 24841, 3455,
                                                                       3485, 12443, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25881, 0, 3,
                                                                       24841, 11879, 24901, 3485,
                                                                       3515, 12503, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25981, 0, 3,
                                                                       24901, 11915, 24961, 3515,
                                                                       3545, 12563, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 26081, 0, 3,
                                                                       24961, 11951, 25021, 3545,
                                                                       3575, 12623, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 26181, 0, 3,
                                                                       25081, 12023, 25181, 3635,
                                                                       3680, 12683, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 26331, 0, 3,
                                                                       25181, 12083, 25281, 3680,
                                                                       3725, 12773, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 26481, 0, 3,
                                                                       25281, 12143, 25381, 3725,
                                                                       3770, 12863, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 26631, 0, 3,
                                                                       25381, 12203, 25481, 3770,
                                                                       3815, 12953, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 26781, 0, 3,
                                                                       25481, 12263, 25581, 3815,
                                                                       3860, 13043, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 26931, 0, 3,
                                                                       25581, 12323, 25681, 3860,
                                                                       3905, 13133, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 27081, 0, 3,
                                                                       25681, 12383, 25781, 3905,
                                                                       3950, 13223, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 27231, 0, 3,
                                                                       25781, 12443, 25881, 3950,
                                                                       3995, 13313, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 27381, 0, 3,
                                                                       25881, 12503, 25981, 3995,
                                                                       4040, 13403, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 27531, 0, 3,
                                                                       25981, 12563, 26081, 4040,
                                                                       4085, 13493, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 27681, 0, 3,
                                                                       26181, 12683, 26331, 4175,
                                                                       4238, 13583, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 27891, 0, 3,
                                                                       26331, 12773, 26481, 4238,
                                                                       4301, 13709, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 28101, 0, 3,
                                                                       26481, 12863, 26631, 4301,
                                                                       4364, 13835, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 28311, 0, 3,
                                                                       26631, 12953, 26781, 4364,
                                                                       4427, 13961, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 28521, 0, 3,
                                                                       26781, 13043, 26931, 4427,
                                                                       4490, 14087, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 28731, 0, 3,
                                                                       26931, 13133, 27081, 4490,
                                                                       4553, 14213, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 28941, 0, 3,
                                                                       27081, 13223, 27231, 4553,
                                                                       4616, 14339, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 29151, 0, 3,
                                                                       27231, 13313, 27381, 4616,
                                                                       4679, 14465, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 29361, 0, 3,
                                                                       27381, 13403, 27531, 4679,
                                                                       4742, 14591, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 29571, 0, 3,
                                                                       27681, 13583, 27891, 4868,
                                                                       4952, 14717, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 29851, 0, 3,
                                                                       27891, 13709, 28101, 4952,
                                                                       5036, 14885, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 30131, 0, 3,
                                                                       28101, 13835, 28311, 5036,
                                                                       5120, 15053, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 30411, 0, 3,
                                                                       28311, 13961, 28521, 5120,
                                                                       5204, 15221, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 30691, 0, 3,
                                                                       28521, 14087, 28731, 5204,
                                                                       5288, 15389, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 30971, 0, 3,
                                                                       28731, 14213, 28941, 5288,
                                                                       5372, 15557, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 31251, 0, 3,
                                                                       28941, 14339, 29151, 5372,
                                                                       5456, 15725, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 31531, 0, 3,
                                                                       29151, 14465, 29361, 5456,
                                                                       5540, 15893, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 31811, 0, 3,
                                                                       29571, 14717, 29851, 5708,
                                                                       5816, 16061, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 32171, 0, 3,
                                                                       29851, 14885, 30131, 5816,
                                                                       5924, 16277, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 32531, 0, 3,
                                                                       30131, 15053, 30411, 5924,
                                                                       6032, 16493, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 32891, 0, 3,
                                                                       30411, 15221, 30691, 6032,
                                                                       6140, 16709, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 33251, 0, 3,
                                                                       30691, 15389, 30971, 6140,
                                                                       6248, 16925, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 33611, 0, 3,
                                                                       30971, 15557, 31251, 6248,
                                                                       6356, 17141, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 33971, 0, 3,
                                                                       31251, 15725, 31531, 6356,
                                                                       6464, 17357, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 34331, 0, 3,
                                                                       31811, 16061, 32171, 6680,
                                                                       6815, 17573, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 34781, 0, 3,
                                                                       32171, 16277, 32531, 6815,
                                                                       6950, 17843, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 35231, 0, 3,
                                                                       32531, 16493, 32891, 6950,
                                                                       7085, 18113, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 35681, 0, 3,
                                                                       32891, 16709, 33251, 7085,
                                                                       7220, 18383, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 36131, 0, 3,
                                                                       33251, 16925, 33611, 7220,
                                                                       7355, 18653, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 36581, 0, 3,
                                                                       33611, 17141, 33971, 7355,
                                                                       7490, 18923, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 37031, 0, 3,
                                                                       34331, 17573, 34781, 7760,
                                                                       7925, 19193, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 37581, 0, 3,
                                                                       34781, 17843, 35231, 7925,
                                                                       8090, 19523, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 38131, 0, 3,
                                                                       35231, 18113, 35681, 8090,
                                                                       8255, 19853, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 38681, 0, 3,
                                                                       35681, 18383, 36131, 8255,
                                                                       8420, 20183, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 39231, 0, 3,
                                                                       36131, 18653, 36581, 8420,
                                                                       8585, 20513, ncols, gamma,
                                                                       p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 39781, 0, 3,
                                                                       37031, 19193, 37581, 8915,
                                                                       9113, 20843, ncols, gamma,
                                                                       p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 40441, 0, 3,
                                                                       37581, 19523, 38131, 9113,
                                                                       9311, 21239, ncols, gamma,
                                                                       p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 41101, 0, 3,
                                                                       38131, 19853, 38681, 9311,
                                                                       9509, 21635, ncols, gamma,
                                                                       p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 41761, 0, 3,
                                                                       38681, 20183, 39231, 9509,
                                                                       9707, 22031, ncols, gamma,
                                                                       p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 42421, 0, 3,
                                                                       39781, 20843, 40441,
                                                                       10103, 10337, 22427,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 43201, 0, 3,
                                                                       40441, 21239, 41101,
                                                                       10337, 10571, 22895,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 43981, 0, 3,
                                                                       41101, 21635, 41761,
                                                                       10571, 10805, 23363,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44761, 3, 11273,
                                                                       11279, 23851, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44776, 3, 11279,
                                                                       11285, 23861, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44791, 3, 11285,
                                                                       11291, 23871, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44806, 3, 11291,
                                                                       11297, 23881, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44821, 3, 11297,
                                                                       11303, 23891, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44836, 3, 11303,
                                                                       11309, 23901, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44851, 3, 11309,
                                                                       11315, 23911, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44866, 3, 11315,
                                                                       11321, 23921, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44881, 3, 11321,
                                                                       11327, 23931, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44896, 3, 11327,
                                                                       11333, 23941, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44911, 3, 11333,
                                                                       11339, 23951, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44926, 3, 11339,
                                                                       11345, 23961, ncols,
                                                                       gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 44941, 0, 3,
                                                                       44761, 23851, 44776,
                                                                       11357, 11375, 24031,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 44986, 0, 3,
                                                                       44776, 23861, 44791,
                                                                       11375, 11393, 24061,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 45031, 0, 3,
                                                                       44791, 23871, 44806,
                                                                       11393, 11411, 24091,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 45076, 0, 3,
                                                                       44806, 23881, 44821,
                                                                       11411, 11429, 24121,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 45121, 0, 3,
                                                                       44821, 23891, 44836,
                                                                       11429, 11447, 24151,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 45166, 0, 3,
                                                                       44836, 23901, 44851,
                                                                       11447, 11465, 24181,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 45211, 0, 3,
                                                                       44851, 23911, 44866,
                                                                       11465, 11483, 24211,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 45256, 0, 3,
                                                                       44866, 23921, 44881,
                                                                       11483, 11501, 24241,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 45301, 0, 3,
                                                                       44881, 23931, 44896,
                                                                       11501, 11519, 24271,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 45346, 0, 3,
                                                                       44896, 23941, 44911,
                                                                       11519, 11537, 24301,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 45391, 0, 3,
                                                                       44911, 23951, 44926,
                                                                       11537, 11555, 24331,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 45436, 0, 3,
                                                                       44941, 24031, 44986,
                                                                       11591, 11627, 24481,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 45526, 0, 3,
                                                                       44986, 24061, 45031,
                                                                       11627, 11663, 24541,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 45616, 0, 3,
                                                                       45031, 24091, 45076,
                                                                       11663, 11699, 24601,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 45706, 0, 3,
                                                                       45076, 24121, 45121,
                                                                       11699, 11735, 24661,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 45796, 0, 3,
                                                                       45121, 24151, 45166,
                                                                       11735, 11771, 24721,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 45886, 0, 3,
                                                                       45166, 24181, 45211,
                                                                       11771, 11807, 24781,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 45976, 0, 3,
                                                                       45211, 24211, 45256,
                                                                       11807, 11843, 24841,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 46066, 0, 3,
                                                                       45256, 24241, 45301,
                                                                       11843, 11879, 24901,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 46156, 0, 3,
                                                                       45301, 24271, 45346,
                                                                       11879, 11915, 24961,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 46246, 0, 3,
                                                                       45346, 24301, 45391,
                                                                       11915, 11951, 25021,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 46336, 0, 3,
                                                                       45436, 24481, 45526,
                                                                       12023, 12083, 25281,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 46486, 0, 3,
                                                                       45526, 24541, 45616,
                                                                       12083, 12143, 25381,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 46636, 0, 3,
                                                                       45616, 24601, 45706,
                                                                       12143, 12203, 25481,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 46786, 0, 3,
                                                                       45706, 24661, 45796,
                                                                       12203, 12263, 25581,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 46936, 0, 3,
                                                                       45796, 24721, 45886,
                                                                       12263, 12323, 25681,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 47086, 0, 3,
                                                                       45886, 24781, 45976,
                                                                       12323, 12383, 25781,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 47236, 0, 3,
                                                                       45976, 24841, 46066,
                                                                       12383, 12443, 25881,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 47386, 0, 3,
                                                                       46066, 24901, 46156,
                                                                       12443, 12503, 25981,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 47536, 0, 3,
                                                                       46156, 24961, 46246,
                                                                       12503, 12563, 26081,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 47686, 0, 3,
                                                                       46336, 25281, 46486,
                                                                       12683, 12773, 26481,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 47911, 0, 3,
                                                                       46486, 25381, 46636,
                                                                       12773, 12863, 26631,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 48136, 0, 3,
                                                                       46636, 25481, 46786,
                                                                       12863, 12953, 26781,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 48361, 0, 3,
                                                                       46786, 25581, 46936,
                                                                       12953, 13043, 26931,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 48586, 0, 3,
                                                                       46936, 25681, 47086,
                                                                       13043, 13133, 27081,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 48811, 0, 3,
                                                                       47086, 25781, 47236,
                                                                       13133, 13223, 27231,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 49036, 0, 3,
                                                                       47236, 25881, 47386,
                                                                       13223, 13313, 27381,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 49261, 0, 3,
                                                                       47386, 25981, 47536,
                                                                       13313, 13403, 27531,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 49486, 0, 3,
                                                                       47686, 26481, 47911,
                                                                       13583, 13709, 28101,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 49801, 0, 3,
                                                                       47911, 26631, 48136,
                                                                       13709, 13835, 28311,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 50116, 0, 3,
                                                                       48136, 26781, 48361,
                                                                       13835, 13961, 28521,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 50431, 0, 3,
                                                                       48361, 26931, 48586,
                                                                       13961, 14087, 28731,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 50746, 0, 3,
                                                                       48586, 27081, 48811,
                                                                       14087, 14213, 28941,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 51061, 0, 3,
                                                                       48811, 27231, 49036,
                                                                       14213, 14339, 29151,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 51376, 0, 3,
                                                                       49036, 27381, 49261,
                                                                       14339, 14465, 29361,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 51691, 0, 3,
                                                                       49486, 28101, 49801,
                                                                       14717, 14885, 30131,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 52111, 0, 3,
                                                                       49801, 28311, 50116,
                                                                       14885, 15053, 30411,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 52531, 0, 3,
                                                                       50116, 28521, 50431,
                                                                       15053, 15221, 30691,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 52951, 0, 3,
                                                                       50431, 28731, 50746,
                                                                       15221, 15389, 30971,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 53371, 0, 3,
                                                                       50746, 28941, 51061,
                                                                       15389, 15557, 31251,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 53791, 0, 3,
                                                                       51061, 29151, 51376,
                                                                       15557, 15725, 31531,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 54211, 0, 3,
                                                                       51691, 30131, 52111,
                                                                       16061, 16277, 32531,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 54751, 0, 3,
                                                                       52111, 30411, 52531,
                                                                       16277, 16493, 32891,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 55291, 0, 3,
                                                                       52531, 30691, 52951,
                                                                       16493, 16709, 33251,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 55831, 0, 3,
                                                                       52951, 30971, 53371,
                                                                       16709, 16925, 33611,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 56371, 0, 3,
                                                                       53371, 31251, 53791,
                                                                       16925, 17141, 33971,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 56911, 0, 3,
                                                                       54211, 32531, 54751,
                                                                       17573, 17843, 35231,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 57586, 0, 3,
                                                                       54751, 32891, 55291,
                                                                       17843, 18113, 35681,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 58261, 0, 3,
                                                                       55291, 33251, 55831,
                                                                       18113, 18383, 36131,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 58936, 0, 3,
                                                                       55831, 33611, 56371,
                                                                       18383, 18653, 36581,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 59611, 0, 3,
                                                                       56911, 35231, 57586,
                                                                       19193, 19523, 38131,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 60436, 0, 3,
                                                                       57586, 35681, 58261,
                                                                       19523, 19853, 38681,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 61261, 0, 3,
                                                                       58261, 36131, 58936,
                                                                       19853, 20183, 39231,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 62086, 0, 3,
                                                                       59611, 38131, 60436,
                                                                       20843, 21239, 41101,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 63076, 0, 3,
                                                                       60436, 38681, 61261,
                                                                       21239, 21635, 41761,
                                                                       ncols, gamma, p, q);

                    compute_prim_sog_three_center_electron_repulsion_0(buffer, 64066, 0, 3,
                                                                       62086, 41101, 63076,
                                                                       22427, 22895, 43981,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65236, 3, 23831,
                                                                       23841, 44761, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65257, 3, 23841,
                                                                       23851, 44776, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65278, 3, 23851,
                                                                       23861, 44791, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65299, 3, 23861,
                                                                       23871, 44806, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65320, 3, 23871,
                                                                       23881, 44821, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65341, 3, 23881,
                                                                       23891, 44836, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65362, 3, 23891,
                                                                       23901, 44851, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65383, 3, 23901,
                                                                       23911, 44866, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65404, 3, 23911,
                                                                       23921, 44881, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65425, 3, 23921,
                                                                       23931, 44896, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65446, 3, 23931,
                                                                       23941, 44911, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65467, 3, 23941,
                                                                       23951, 44926, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 65488, 0, 3,
                                                                       65236, 44761, 65257,
                                                                       23971, 24001, 44941,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 65551, 0, 3,
                                                                       65257, 44776, 65278,
                                                                       24001, 24031, 44986,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 65614, 0, 3,
                                                                       65278, 44791, 65299,
                                                                       24031, 24061, 45031,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 65677, 0, 3,
                                                                       65299, 44806, 65320,
                                                                       24061, 24091, 45076,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 65740, 0, 3,
                                                                       65320, 44821, 65341,
                                                                       24091, 24121, 45121,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 65803, 0, 3,
                                                                       65341, 44836, 65362,
                                                                       24121, 24151, 45166,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 65866, 0, 3,
                                                                       65362, 44851, 65383,
                                                                       24151, 24181, 45211,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 65929, 0, 3,
                                                                       65383, 44866, 65404,
                                                                       24181, 24211, 45256,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 65992, 0, 3,
                                                                       65404, 44881, 65425,
                                                                       24211, 24241, 45301,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 66055, 0, 3,
                                                                       65425, 44896, 65446,
                                                                       24241, 24271, 45346,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 66118, 0, 3,
                                                                       65446, 44911, 65467,
                                                                       24271, 24301, 45391,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 66181, 0, 3,
                                                                       65488, 44941, 65551,
                                                                       24361, 24421, 45436,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 66307, 0, 3,
                                                                       65551, 44986, 65614,
                                                                       24421, 24481, 45526,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 66433, 0, 3,
                                                                       65614, 45031, 65677,
                                                                       24481, 24541, 45616,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 66559, 0, 3,
                                                                       65677, 45076, 65740,
                                                                       24541, 24601, 45706,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 66685, 0, 3,
                                                                       65740, 45121, 65803,
                                                                       24601, 24661, 45796,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 66811, 0, 3,
                                                                       65803, 45166, 65866,
                                                                       24661, 24721, 45886,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 66937, 0, 3,
                                                                       65866, 45211, 65929,
                                                                       24721, 24781, 45976,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 67063, 0, 3,
                                                                       65929, 45256, 65992,
                                                                       24781, 24841, 46066,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 67189, 0, 3,
                                                                       65992, 45301, 66055,
                                                                       24841, 24901, 46156,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 67315, 0, 3,
                                                                       66055, 45346, 66118,
                                                                       24901, 24961, 46246,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 67441, 0, 3,
                                                                       66181, 45436, 66307,
                                                                       25081, 25181, 46336,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 67651, 0, 3,
                                                                       66307, 45526, 66433,
                                                                       25181, 25281, 46486,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 67861, 0, 3,
                                                                       66433, 45616, 66559,
                                                                       25281, 25381, 46636,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 68071, 0, 3,
                                                                       66559, 45706, 66685,
                                                                       25381, 25481, 46786,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 68281, 0, 3,
                                                                       66685, 45796, 66811,
                                                                       25481, 25581, 46936,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 68491, 0, 3,
                                                                       66811, 45886, 66937,
                                                                       25581, 25681, 47086,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 68701, 0, 3,
                                                                       66937, 45976, 67063,
                                                                       25681, 25781, 47236,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 68911, 0, 3,
                                                                       67063, 46066, 67189,
                                                                       25781, 25881, 47386,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 69121, 0, 3,
                                                                       67189, 46156, 67315,
                                                                       25881, 25981, 47536,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 69331, 0, 3,
                                                                       67441, 46336, 67651,
                                                                       26181, 26331, 47686,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 69646, 0, 3,
                                                                       67651, 46486, 67861,
                                                                       26331, 26481, 47911,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 69961, 0, 3,
                                                                       67861, 46636, 68071,
                                                                       26481, 26631, 48136,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 70276, 0, 3,
                                                                       68071, 46786, 68281,
                                                                       26631, 26781, 48361,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 70591, 0, 3,
                                                                       68281, 46936, 68491,
                                                                       26781, 26931, 48586,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 70906, 0, 3,
                                                                       68491, 47086, 68701,
                                                                       26931, 27081, 48811,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 71221, 0, 3,
                                                                       68701, 47236, 68911,
                                                                       27081, 27231, 49036,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 71536, 0, 3,
                                                                       68911, 47386, 69121,
                                                                       27231, 27381, 49261,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 71851, 0, 3,
                                                                       69331, 47686, 69646,
                                                                       27681, 27891, 49486,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 72292, 0, 3,
                                                                       69646, 47911, 69961,
                                                                       27891, 28101, 49801,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 72733, 0, 3,
                                                                       69961, 48136, 70276,
                                                                       28101, 28311, 50116,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 73174, 0, 3,
                                                                       70276, 48361, 70591,
                                                                       28311, 28521, 50431,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 73615, 0, 3,
                                                                       70591, 48586, 70906,
                                                                       28521, 28731, 50746,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 74056, 0, 3,
                                                                       70906, 48811, 71221,
                                                                       28731, 28941, 51061,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 74497, 0, 3,
                                                                       71221, 49036, 71536,
                                                                       28941, 29151, 51376,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 74938, 0, 3,
                                                                       71851, 49486, 72292,
                                                                       29571, 29851, 51691,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 75526, 0, 3,
                                                                       72292, 49801, 72733,
                                                                       29851, 30131, 52111,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 76114, 0, 3,
                                                                       72733, 50116, 73174,
                                                                       30131, 30411, 52531,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 76702, 0, 3,
                                                                       73174, 50431, 73615,
                                                                       30411, 30691, 52951,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 77290, 0, 3,
                                                                       73615, 50746, 74056,
                                                                       30691, 30971, 53371,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 77878, 0, 3,
                                                                       74056, 51061, 74497,
                                                                       30971, 31251, 53791,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 78466, 0, 3,
                                                                       74938, 51691, 75526,
                                                                       31811, 32171, 54211,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 79222, 0, 3,
                                                                       75526, 52111, 76114,
                                                                       32171, 32531, 54751,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 79978, 0, 3,
                                                                       76114, 52531, 76702,
                                                                       32531, 32891, 55291,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 80734, 0, 3,
                                                                       76702, 52951, 77290,
                                                                       32891, 33251, 55831,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 81490, 0, 3,
                                                                       77290, 53371, 77878,
                                                                       33251, 33611, 56371,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 82246, 0, 3,
                                                                       78466, 54211, 79222,
                                                                       34331, 34781, 56911,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 83191, 0, 3,
                                                                       79222, 54751, 79978,
                                                                       34781, 35231, 57586,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 84136, 0, 3,
                                                                       79978, 55291, 80734,
                                                                       35231, 35681, 58261,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 85081, 0, 3,
                                                                       80734, 55831, 81490,
                                                                       35681, 36131, 58936,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 86026, 0, 3,
                                                                       82246, 56911, 83191,
                                                                       37031, 37581, 59611,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 87181, 0, 3,
                                                                       83191, 57586, 84136,
                                                                       37581, 38131, 60436,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 88336, 0, 3,
                                                                       84136, 58261, 85081,
                                                                       38131, 38681, 61261,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 89491, 0, 3,
                                                                       86026, 59611, 87181,
                                                                       39781, 40441, 62086,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 90877, 0, 3,
                                                                       87181, 60436, 88336,
                                                                       40441, 41101, 63076,
                                                                       ncols, gamma, p, q);

                    compute_prim_soh_three_center_electron_repulsion_0(buffer, 92263, 0, 3,
                                                                       89491, 62086, 90877,
                                                                       42421, 43201, 64066,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 93901, 74938, 588, ncols);

                    simdfunc::contract_primitives(buffer, 94797, 78466, 756, ncols);

                    simdfunc::contract_primitives(buffer, 95949, 82246, 945, ncols);

                    simdfunc::contract_primitives(buffer, 97389, 86026, 1155, ncols);

                    simdfunc::contract_primitives(buffer, 99149, 89491, 1386, ncols);

                    simdfunc::contract_primitives(buffer, 101261, 92263, 1638, ncols);
                }
            }
        }

        simdtrf::transform_h_inner(buffer, 94489, 93901, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 95553, 94797, 36, 1, nmax);

        simdtrf::transform_h_inner(buffer, 96894, 95949, 45, 1, nmax);

        simdtrf::transform_h_inner(buffer, 98544, 97389, 55, 1, nmax);

        simdtrf::transform_h_inner(buffer, 100535, 99149, 66, 1, nmax);

        simdtrf::transform_h_inner(buffer, 102899, 101261, 78, 1, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 103757, 94489, 95553, 11, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 104681, 95553, 96894, 11, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 105869, 96894, 98544, 11, nmax);

        simdtrf::compute_hrr_pm(buffer, coordinates, 107354, 98544, 100535, 11, nmax);

        simdtrf::compute_hrr_pn(buffer, coordinates, 109169, 100535, 102899, 11, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 111347, 103757, 104681, 11, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 113195, 104681, 105869, 11, nmax);

        simdtrf::compute_hrr_dl(buffer, coordinates, 115571, 105869, 107354, 11, nmax);

        simdtrf::compute_hrr_dm(buffer, coordinates, 118541, 107354, 109169, 11, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 122171, 111347, 113195, 11, nmax);

        simdtrf::compute_hrr_fk(buffer, coordinates, 125251, 113195, 115571, 11, nmax);

        simdtrf::compute_hrr_fl(buffer, coordinates, 129211, 115571, 118541, 11, nmax);

        simdtrf::compute_hrr_gi(buffer, coordinates, 134161, 122171, 125251, 11, nmax);

        simdtrf::compute_hrr_gk(buffer, coordinates, 138781, 125251, 129211, 11, nmax);

        simdtrf::compute_hrr_hi(buffer, coordinates, 144721, 134161, 138781, 11, nmax);

        simdtrf::transform_i_inner(buffer, 151189, 144721, 21, 11, nmax);

        simdtrf::transform_h_outer(values + n * npairs, nvalues, buffer, 151189, 143, nmax);
    }

    for (size_t m = 0; m < 1573; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
