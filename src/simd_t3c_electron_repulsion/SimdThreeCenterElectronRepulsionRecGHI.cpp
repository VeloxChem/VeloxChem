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


#include "SimdThreeCenterElectronRepulsionRecGHI.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSID.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSII.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSI.hpp"
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
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_ghi_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_ghi_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 112202, 0, 0, dimensions);

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
        simdfunc::prepare_buffer(buffer, 112202, 80307, 6870, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 6, 3, 15,
                                                             ncols, fj, mu, fq);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2042, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2045, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2048, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2051, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2054, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2057, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2060, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2063, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2066, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2069, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2072, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2075, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2078, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2081, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2084, 3, 9, 29,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2093, 3, 10, 32,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2102, 3, 11, 35,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2111, 3, 12, 38,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2120, 3, 13, 41,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2129, 3, 14, 44,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2138, 3, 15, 47,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2147, 3, 16, 50,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2156, 3, 17, 53,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2165, 3, 18, 56,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2174, 3, 19, 59,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2183, 3, 20, 62,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2192, 3, 21, 65,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2201, 3, 29, 80,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2219, 3, 32, 86,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2237, 3, 35, 92,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2255, 3, 38, 98,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2273, 3, 41, 104,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2291, 3, 44, 110,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2309, 3, 47, 116,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2327, 3, 50, 122,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2345, 3, 53, 128,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2363, 3, 56, 134,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2381, 3, 59, 140,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2399, 3, 62, 146,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2417, 3, 80, 172,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2447, 3, 86, 182,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2477, 3, 92, 192,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2507, 3, 98, 202,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2537, 3, 104, 212,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2567, 3, 110, 222,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2597, 3, 116, 232,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2627, 3, 122, 242,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2657, 3, 128, 252,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2687, 3, 134, 262,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2717, 3, 140, 272,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2747, 3, 172, 312,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2792, 3, 182, 327,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2837, 3, 192, 342,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2882, 3, 202, 357,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2927, 3, 212, 372,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2972, 3, 222, 387,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3017, 3, 232, 402,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3062, 3, 242, 417,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3107, 3, 252, 432,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3152, 3, 262, 447,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3197, 3, 312, 504,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3260, 3, 327, 525,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3323, 3, 342, 546,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3386, 3, 357, 567,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3449, 3, 372, 588,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3512, 3, 387, 609,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3575, 3, 402, 630,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3638, 3, 417, 651,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3701, 3, 432, 672,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3764, 3, 504, 749,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3848, 3, 525, 777,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3932, 3, 546, 805,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4016, 3, 567, 833,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4100, 3, 588, 861,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4184, 3, 609, 889,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4268, 3, 630, 917,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4352, 3, 651, 945,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4436, 3, 749,
                                                                       1045, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4544, 3, 777,
                                                                       1081, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4652, 3, 805,
                                                                       1117, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4760, 3, 833,
                                                                       1153, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4868, 3, 861,
                                                                       1189, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4976, 3, 889,
                                                                       1225, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5084, 3, 917,
                                                                       1261, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5192, 3, 1045,
                                                                       1387, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5327, 3, 1081,
                                                                       1432, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5462, 3, 1117,
                                                                       1477, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5597, 3, 1153,
                                                                       1522, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5732, 3, 1189,
                                                                       1567, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5867, 3, 1225,
                                                                       1612, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 6002, 3, 1387,
                                                                       1767, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 6167, 3, 1432,
                                                                       1822, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 6332, 3, 1477,
                                                                       1877, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 6497, 3, 1522,
                                                                       1932, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 6662, 3, 1567,
                                                                       1987, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6827, 3, 7, 8,
                                                                       2042, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6833, 3, 8, 9,
                                                                       2045, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6839, 3, 9, 10,
                                                                       2048, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6845, 3, 10, 11,
                                                                       2051, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6851, 3, 11, 12,
                                                                       2054, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6857, 3, 12, 13,
                                                                       2057, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6863, 3, 13, 14,
                                                                       2060, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6869, 3, 14, 15,
                                                                       2063, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6875, 3, 15, 16,
                                                                       2066, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6881, 3, 16, 17,
                                                                       2069, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6887, 3, 17, 18,
                                                                       2072, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6893, 3, 18, 19,
                                                                       2075, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6899, 3, 19, 20,
                                                                       2078, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6905, 3, 20, 21,
                                                                       2081, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6911, 0, 3, 6827,
                                                                       2042, 6833, 2084, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6929, 0, 3, 6833,
                                                                       2045, 6839, 2093, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6947, 0, 3, 6839,
                                                                       2048, 6845, 2102, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6965, 0, 3, 6845,
                                                                       2051, 6851, 2111, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6983, 0, 3, 6851,
                                                                       2054, 6857, 2120, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7001, 0, 3, 6857,
                                                                       2057, 6863, 2129, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7019, 0, 3, 6863,
                                                                       2060, 6869, 2138, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7037, 0, 3, 6869,
                                                                       2063, 6875, 2147, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7055, 0, 3, 6875,
                                                                       2066, 6881, 2156, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7073, 0, 3, 6881,
                                                                       2069, 6887, 2165, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7091, 0, 3, 6887,
                                                                       2072, 6893, 2174, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7109, 0, 3, 6893,
                                                                       2075, 6899, 2183, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7127, 0, 3, 6899,
                                                                       2078, 6905, 2192, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7145, 0, 3, 6911,
                                                                       2084, 6929, 68, 74, 2201,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7181, 0, 3, 6929,
                                                                       2093, 6947, 74, 80, 2219,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7217, 0, 3, 6947,
                                                                       2102, 6965, 80, 86, 2237,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7253, 0, 3, 6965,
                                                                       2111, 6983, 86, 92, 2255,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7289, 0, 3, 6983,
                                                                       2120, 7001, 92, 98, 2273,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7325, 0, 3, 7001,
                                                                       2129, 7019, 98, 104, 2291,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7361, 0, 3, 7019,
                                                                       2138, 7037, 104, 110,
                                                                       2309, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7397, 0, 3, 7037,
                                                                       2147, 7055, 110, 116,
                                                                       2327, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7433, 0, 3, 7055,
                                                                       2156, 7073, 116, 122,
                                                                       2345, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7469, 0, 3, 7073,
                                                                       2165, 7091, 122, 128,
                                                                       2363, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7505, 0, 3, 7091,
                                                                       2174, 7109, 128, 134,
                                                                       2381, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7541, 0, 3, 7109,
                                                                       2183, 7127, 134, 140,
                                                                       2399, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7577, 0, 3, 7145,
                                                                       2201, 7181, 152, 162,
                                                                       2417, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7637, 0, 3, 7181,
                                                                       2219, 7217, 162, 172,
                                                                       2447, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7697, 0, 3, 7217,
                                                                       2237, 7253, 172, 182,
                                                                       2477, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7757, 0, 3, 7253,
                                                                       2255, 7289, 182, 192,
                                                                       2507, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7817, 0, 3, 7289,
                                                                       2273, 7325, 192, 202,
                                                                       2537, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7877, 0, 3, 7325,
                                                                       2291, 7361, 202, 212,
                                                                       2567, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7937, 0, 3, 7361,
                                                                       2309, 7397, 212, 222,
                                                                       2597, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7997, 0, 3, 7397,
                                                                       2327, 7433, 222, 232,
                                                                       2627, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8057, 0, 3, 7433,
                                                                       2345, 7469, 232, 242,
                                                                       2657, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8117, 0, 3, 7469,
                                                                       2363, 7505, 242, 252,
                                                                       2687, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8177, 0, 3, 7505,
                                                                       2381, 7541, 252, 262,
                                                                       2717, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8237, 0, 3, 7577,
                                                                       2417, 7637, 282, 297,
                                                                       2747, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8327, 0, 3, 7637,
                                                                       2447, 7697, 297, 312,
                                                                       2792, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8417, 0, 3, 7697,
                                                                       2477, 7757, 312, 327,
                                                                       2837, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8507, 0, 3, 7757,
                                                                       2507, 7817, 327, 342,
                                                                       2882, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8597, 0, 3, 7817,
                                                                       2537, 7877, 342, 357,
                                                                       2927, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8687, 0, 3, 7877,
                                                                       2567, 7937, 357, 372,
                                                                       2972, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8777, 0, 3, 7937,
                                                                       2597, 7997, 372, 387,
                                                                       3017, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8867, 0, 3, 7997,
                                                                       2627, 8057, 387, 402,
                                                                       3062, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8957, 0, 3, 8057,
                                                                       2657, 8117, 402, 417,
                                                                       3107, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9047, 0, 3, 8117,
                                                                       2687, 8177, 417, 432,
                                                                       3152, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9137, 0, 3, 8237,
                                                                       2747, 8327, 462, 483,
                                                                       3197, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9263, 0, 3, 8327,
                                                                       2792, 8417, 483, 504,
                                                                       3260, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9389, 0, 3, 8417,
                                                                       2837, 8507, 504, 525,
                                                                       3323, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9515, 0, 3, 8507,
                                                                       2882, 8597, 525, 546,
                                                                       3386, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9641, 0, 3, 8597,
                                                                       2927, 8687, 546, 567,
                                                                       3449, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9767, 0, 3, 8687,
                                                                       2972, 8777, 567, 588,
                                                                       3512, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9893, 0, 3, 8777,
                                                                       3017, 8867, 588, 609,
                                                                       3575, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10019, 0, 3, 8867,
                                                                       3062, 8957, 609, 630,
                                                                       3638, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10145, 0, 3, 8957,
                                                                       3107, 9047, 630, 651,
                                                                       3701, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 10271, 0, 3, 9137,
                                                                       3197, 9263, 693, 721,
                                                                       3764, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 10439, 0, 3, 9263,
                                                                       3260, 9389, 721, 749,
                                                                       3848, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 10607, 0, 3, 9389,
                                                                       3323, 9515, 749, 777,
                                                                       3932, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 10775, 0, 3, 9515,
                                                                       3386, 9641, 777, 805,
                                                                       4016, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 10943, 0, 3, 9641,
                                                                       3449, 9767, 805, 833,
                                                                       4100, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11111, 0, 3, 9767,
                                                                       3512, 9893, 833, 861,
                                                                       4184, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11279, 0, 3, 9893,
                                                                       3575, 10019, 861, 889,
                                                                       4268, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11447, 0, 3,
                                                                       10019, 3638, 10145, 889,
                                                                       917, 4352, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 11615, 0, 3,
                                                                       10271, 3764, 10439, 973,
                                                                       1009, 4436, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 11831, 0, 3,
                                                                       10439, 3848, 10607, 1009,
                                                                       1045, 4544, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 12047, 0, 3,
                                                                       10607, 3932, 10775, 1045,
                                                                       1081, 4652, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 12263, 0, 3,
                                                                       10775, 4016, 10943, 1081,
                                                                       1117, 4760, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 12479, 0, 3,
                                                                       10943, 4100, 11111, 1117,
                                                                       1153, 4868, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 12695, 0, 3,
                                                                       11111, 4184, 11279, 1153,
                                                                       1189, 4976, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 12911, 0, 3,
                                                                       11279, 4268, 11447, 1189,
                                                                       1225, 5084, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 13127, 0, 3,
                                                                       11615, 4436, 11831, 1297,
                                                                       1342, 5192, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 13397, 0, 3,
                                                                       11831, 4544, 12047, 1342,
                                                                       1387, 5327, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 13667, 0, 3,
                                                                       12047, 4652, 12263, 1387,
                                                                       1432, 5462, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 13937, 0, 3,
                                                                       12263, 4760, 12479, 1432,
                                                                       1477, 5597, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 14207, 0, 3,
                                                                       12479, 4868, 12695, 1477,
                                                                       1522, 5732, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 14477, 0, 3,
                                                                       12695, 4976, 12911, 1522,
                                                                       1567, 5867, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 14747, 0, 3,
                                                                       13127, 5192, 13397, 1657,
                                                                       1712, 6002, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 15077, 0, 3,
                                                                       13397, 5327, 13667, 1712,
                                                                       1767, 6167, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 15407, 0, 3,
                                                                       13667, 5462, 13937, 1767,
                                                                       1822, 6332, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 15737, 0, 3,
                                                                       13937, 5597, 14207, 1822,
                                                                       1877, 6497, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 16067, 0, 3,
                                                                       14207, 5732, 14477, 1877,
                                                                       1932, 6662, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16397, 3, 2042,
                                                                       2045, 6839, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16407, 3, 2045,
                                                                       2048, 6845, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16417, 3, 2048,
                                                                       2051, 6851, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16427, 3, 2051,
                                                                       2054, 6857, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16437, 3, 2054,
                                                                       2057, 6863, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16447, 3, 2057,
                                                                       2060, 6869, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16457, 3, 2060,
                                                                       2063, 6875, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16467, 3, 2063,
                                                                       2066, 6881, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16477, 3, 2066,
                                                                       2069, 6887, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16487, 3, 2069,
                                                                       2072, 6893, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16497, 3, 2072,
                                                                       2075, 6899, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16507, 3, 2075,
                                                                       2078, 6905, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16517, 0, 3,
                                                                       16397, 6839, 16407, 6947,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16547, 0, 3,
                                                                       16407, 6845, 16417, 6965,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16577, 0, 3,
                                                                       16417, 6851, 16427, 6983,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16607, 0, 3,
                                                                       16427, 6857, 16437, 7001,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16637, 0, 3,
                                                                       16437, 6863, 16447, 7019,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16667, 0, 3,
                                                                       16447, 6869, 16457, 7037,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16697, 0, 3,
                                                                       16457, 6875, 16467, 7055,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16727, 0, 3,
                                                                       16467, 6881, 16477, 7073,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16757, 0, 3,
                                                                       16477, 6887, 16487, 7091,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16787, 0, 3,
                                                                       16487, 6893, 16497, 7109,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16817, 0, 3,
                                                                       16497, 6899, 16507, 7127,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 16847, 0, 3,
                                                                       16517, 6947, 16547, 2201,
                                                                       2219, 7217, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 16907, 0, 3,
                                                                       16547, 6965, 16577, 2219,
                                                                       2237, 7253, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 16967, 0, 3,
                                                                       16577, 6983, 16607, 2237,
                                                                       2255, 7289, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 17027, 0, 3,
                                                                       16607, 7001, 16637, 2255,
                                                                       2273, 7325, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 17087, 0, 3,
                                                                       16637, 7019, 16667, 2273,
                                                                       2291, 7361, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 17147, 0, 3,
                                                                       16667, 7037, 16697, 2291,
                                                                       2309, 7397, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 17207, 0, 3,
                                                                       16697, 7055, 16727, 2309,
                                                                       2327, 7433, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 17267, 0, 3,
                                                                       16727, 7073, 16757, 2327,
                                                                       2345, 7469, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 17327, 0, 3,
                                                                       16757, 7091, 16787, 2345,
                                                                       2363, 7505, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 17387, 0, 3,
                                                                       16787, 7109, 16817, 2363,
                                                                       2381, 7541, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 17447, 0, 3,
                                                                       16847, 7217, 16907, 2417,
                                                                       2447, 7697, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 17547, 0, 3,
                                                                       16907, 7253, 16967, 2447,
                                                                       2477, 7757, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 17647, 0, 3,
                                                                       16967, 7289, 17027, 2477,
                                                                       2507, 7817, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 17747, 0, 3,
                                                                       17027, 7325, 17087, 2507,
                                                                       2537, 7877, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 17847, 0, 3,
                                                                       17087, 7361, 17147, 2537,
                                                                       2567, 7937, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 17947, 0, 3,
                                                                       17147, 7397, 17207, 2567,
                                                                       2597, 7997, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 18047, 0, 3,
                                                                       17207, 7433, 17267, 2597,
                                                                       2627, 8057, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 18147, 0, 3,
                                                                       17267, 7469, 17327, 2627,
                                                                       2657, 8117, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 18247, 0, 3,
                                                                       17327, 7505, 17387, 2657,
                                                                       2687, 8177, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 18347, 0, 3,
                                                                       17447, 7697, 17547, 2747,
                                                                       2792, 8417, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 18497, 0, 3,
                                                                       17547, 7757, 17647, 2792,
                                                                       2837, 8507, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 18647, 0, 3,
                                                                       17647, 7817, 17747, 2837,
                                                                       2882, 8597, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 18797, 0, 3,
                                                                       17747, 7877, 17847, 2882,
                                                                       2927, 8687, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 18947, 0, 3,
                                                                       17847, 7937, 17947, 2927,
                                                                       2972, 8777, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 19097, 0, 3,
                                                                       17947, 7997, 18047, 2972,
                                                                       3017, 8867, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 19247, 0, 3,
                                                                       18047, 8057, 18147, 3017,
                                                                       3062, 8957, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 19397, 0, 3,
                                                                       18147, 8117, 18247, 3062,
                                                                       3107, 9047, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 19547, 0, 3,
                                                                       18347, 8417, 18497, 3197,
                                                                       3260, 9389, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 19757, 0, 3,
                                                                       18497, 8507, 18647, 3260,
                                                                       3323, 9515, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 19967, 0, 3,
                                                                       18647, 8597, 18797, 3323,
                                                                       3386, 9641, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 20177, 0, 3,
                                                                       18797, 8687, 18947, 3386,
                                                                       3449, 9767, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 20387, 0, 3,
                                                                       18947, 8777, 19097, 3449,
                                                                       3512, 9893, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 20597, 0, 3,
                                                                       19097, 8867, 19247, 3512,
                                                                       3575, 10019, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 20807, 0, 3,
                                                                       19247, 8957, 19397, 3575,
                                                                       3638, 10145, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 21017, 0, 3,
                                                                       19547, 9389, 19757, 3764,
                                                                       3848, 10607, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 21297, 0, 3,
                                                                       19757, 9515, 19967, 3848,
                                                                       3932, 10775, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 21577, 0, 3,
                                                                       19967, 9641, 20177, 3932,
                                                                       4016, 10943, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 21857, 0, 3,
                                                                       20177, 9767, 20387, 4016,
                                                                       4100, 11111, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 22137, 0, 3,
                                                                       20387, 9893, 20597, 4100,
                                                                       4184, 11279, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 22417, 0, 3,
                                                                       20597, 10019, 20807, 4184,
                                                                       4268, 11447, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 22697, 0, 3,
                                                                       21017, 10607, 21297, 4436,
                                                                       4544, 12047, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 23057, 0, 3,
                                                                       21297, 10775, 21577, 4544,
                                                                       4652, 12263, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 23417, 0, 3,
                                                                       21577, 10943, 21857, 4652,
                                                                       4760, 12479, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 23777, 0, 3,
                                                                       21857, 11111, 22137, 4760,
                                                                       4868, 12695, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 24137, 0, 3,
                                                                       22137, 11279, 22417, 4868,
                                                                       4976, 12911, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 24497, 0, 3,
                                                                       22697, 12047, 23057, 5192,
                                                                       5327, 13667, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 24947, 0, 3,
                                                                       23057, 12263, 23417, 5327,
                                                                       5462, 13937, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 25397, 0, 3,
                                                                       23417, 12479, 23777, 5462,
                                                                       5597, 14207, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 25847, 0, 3,
                                                                       23777, 12695, 24137, 5597,
                                                                       5732, 14477, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 26297, 0, 3,
                                                                       24497, 13667, 24947, 6002,
                                                                       6167, 15407, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 26847, 0, 3,
                                                                       24947, 13937, 25397, 6167,
                                                                       6332, 15737, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 27397, 0, 3,
                                                                       25397, 14207, 25847, 6332,
                                                                       6497, 16067, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 27947, 3, 6827,
                                                                       6833, 16397, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 27962, 3, 6833,
                                                                       6839, 16407, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 27977, 3, 6839,
                                                                       6845, 16417, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 27992, 3, 6845,
                                                                       6851, 16427, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 28007, 3, 6851,
                                                                       6857, 16437, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 28022, 3, 6857,
                                                                       6863, 16447, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 28037, 3, 6863,
                                                                       6869, 16457, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 28052, 3, 6869,
                                                                       6875, 16467, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 28067, 3, 6875,
                                                                       6881, 16477, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 28082, 3, 6881,
                                                                       6887, 16487, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 28097, 3, 6887,
                                                                       6893, 16497, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 28112, 3, 6893,
                                                                       6899, 16507, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 28127, 0, 3,
                                                                       27947, 16397, 27962, 6911,
                                                                       6929, 16517, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 28172, 0, 3,
                                                                       27962, 16407, 27977, 6929,
                                                                       6947, 16547, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 28217, 0, 3,
                                                                       27977, 16417, 27992, 6947,
                                                                       6965, 16577, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 28262, 0, 3,
                                                                       27992, 16427, 28007, 6965,
                                                                       6983, 16607, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 28307, 0, 3,
                                                                       28007, 16437, 28022, 6983,
                                                                       7001, 16637, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 28352, 0, 3,
                                                                       28022, 16447, 28037, 7001,
                                                                       7019, 16667, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 28397, 0, 3,
                                                                       28037, 16457, 28052, 7019,
                                                                       7037, 16697, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 28442, 0, 3,
                                                                       28052, 16467, 28067, 7037,
                                                                       7055, 16727, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 28487, 0, 3,
                                                                       28067, 16477, 28082, 7055,
                                                                       7073, 16757, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 28532, 0, 3,
                                                                       28082, 16487, 28097, 7073,
                                                                       7091, 16787, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 28577, 0, 3,
                                                                       28097, 16497, 28112, 7091,
                                                                       7109, 16817, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 28622, 0, 3,
                                                                       28127, 16517, 28172, 7145,
                                                                       7181, 16847, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 28712, 0, 3,
                                                                       28172, 16547, 28217, 7181,
                                                                       7217, 16907, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 28802, 0, 3,
                                                                       28217, 16577, 28262, 7217,
                                                                       7253, 16967, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 28892, 0, 3,
                                                                       28262, 16607, 28307, 7253,
                                                                       7289, 17027, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 28982, 0, 3,
                                                                       28307, 16637, 28352, 7289,
                                                                       7325, 17087, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 29072, 0, 3,
                                                                       28352, 16667, 28397, 7325,
                                                                       7361, 17147, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 29162, 0, 3,
                                                                       28397, 16697, 28442, 7361,
                                                                       7397, 17207, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 29252, 0, 3,
                                                                       28442, 16727, 28487, 7397,
                                                                       7433, 17267, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 29342, 0, 3,
                                                                       28487, 16757, 28532, 7433,
                                                                       7469, 17327, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 29432, 0, 3,
                                                                       28532, 16787, 28577, 7469,
                                                                       7505, 17387, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 29522, 0, 3,
                                                                       28622, 16847, 28712, 7577,
                                                                       7637, 17447, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 29672, 0, 3,
                                                                       28712, 16907, 28802, 7637,
                                                                       7697, 17547, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 29822, 0, 3,
                                                                       28802, 16967, 28892, 7697,
                                                                       7757, 17647, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 29972, 0, 3,
                                                                       28892, 17027, 28982, 7757,
                                                                       7817, 17747, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 30122, 0, 3,
                                                                       28982, 17087, 29072, 7817,
                                                                       7877, 17847, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 30272, 0, 3,
                                                                       29072, 17147, 29162, 7877,
                                                                       7937, 17947, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 30422, 0, 3,
                                                                       29162, 17207, 29252, 7937,
                                                                       7997, 18047, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 30572, 0, 3,
                                                                       29252, 17267, 29342, 7997,
                                                                       8057, 18147, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 30722, 0, 3,
                                                                       29342, 17327, 29432, 8057,
                                                                       8117, 18247, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 30872, 0, 3,
                                                                       29522, 17447, 29672, 8237,
                                                                       8327, 18347, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 31097, 0, 3,
                                                                       29672, 17547, 29822, 8327,
                                                                       8417, 18497, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 31322, 0, 3,
                                                                       29822, 17647, 29972, 8417,
                                                                       8507, 18647, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 31547, 0, 3,
                                                                       29972, 17747, 30122, 8507,
                                                                       8597, 18797, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 31772, 0, 3,
                                                                       30122, 17847, 30272, 8597,
                                                                       8687, 18947, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 31997, 0, 3,
                                                                       30272, 17947, 30422, 8687,
                                                                       8777, 19097, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 32222, 0, 3,
                                                                       30422, 18047, 30572, 8777,
                                                                       8867, 19247, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 32447, 0, 3,
                                                                       30572, 18147, 30722, 8867,
                                                                       8957, 19397, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 32672, 0, 3,
                                                                       30872, 18347, 31097, 9137,
                                                                       9263, 19547, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 32987, 0, 3,
                                                                       31097, 18497, 31322, 9263,
                                                                       9389, 19757, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 33302, 0, 3,
                                                                       31322, 18647, 31547, 9389,
                                                                       9515, 19967, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 33617, 0, 3,
                                                                       31547, 18797, 31772, 9515,
                                                                       9641, 20177, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 33932, 0, 3,
                                                                       31772, 18947, 31997, 9641,
                                                                       9767, 20387, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 34247, 0, 3,
                                                                       31997, 19097, 32222, 9767,
                                                                       9893, 20597, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 34562, 0, 3,
                                                                       32222, 19247, 32447, 9893,
                                                                       10019, 20807, ncols,
                                                                       gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 34877, 0, 3,
                                                                       32672, 19547, 32987,
                                                                       10271, 10439, 21017,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 35297, 0, 3,
                                                                       32987, 19757, 33302,
                                                                       10439, 10607, 21297,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 35717, 0, 3,
                                                                       33302, 19967, 33617,
                                                                       10607, 10775, 21577,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 36137, 0, 3,
                                                                       33617, 20177, 33932,
                                                                       10775, 10943, 21857,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 36557, 0, 3,
                                                                       33932, 20387, 34247,
                                                                       10943, 11111, 22137,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 36977, 0, 3,
                                                                       34247, 20597, 34562,
                                                                       11111, 11279, 22417,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 37397, 0, 3,
                                                                       34877, 21017, 35297,
                                                                       11615, 11831, 22697,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 37937, 0, 3,
                                                                       35297, 21297, 35717,
                                                                       11831, 12047, 23057,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 38477, 0, 3,
                                                                       35717, 21577, 36137,
                                                                       12047, 12263, 23417,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 39017, 0, 3,
                                                                       36137, 21857, 36557,
                                                                       12263, 12479, 23777,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 39557, 0, 3,
                                                                       36557, 22137, 36977,
                                                                       12479, 12695, 24137,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 40097, 0, 3,
                                                                       37397, 22697, 37937,
                                                                       13127, 13397, 24497,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 40772, 0, 3,
                                                                       37937, 23057, 38477,
                                                                       13397, 13667, 24947,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 41447, 0, 3,
                                                                       38477, 23417, 39017,
                                                                       13667, 13937, 25397,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 42122, 0, 3,
                                                                       39017, 23777, 39557,
                                                                       13937, 14207, 25847,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 42797, 0, 3,
                                                                       40097, 24497, 40772,
                                                                       14747, 15077, 26297,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 43622, 0, 3,
                                                                       40772, 24947, 41447,
                                                                       15077, 15407, 26847,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 44447, 0, 3,
                                                                       41447, 25397, 42122,
                                                                       15407, 15737, 27397,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 45272, 3, 16397,
                                                                       16407, 27977, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 45293, 3, 16407,
                                                                       16417, 27992, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 45314, 3, 16417,
                                                                       16427, 28007, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 45335, 3, 16427,
                                                                       16437, 28022, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 45356, 3, 16437,
                                                                       16447, 28037, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 45377, 3, 16447,
                                                                       16457, 28052, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 45398, 3, 16457,
                                                                       16467, 28067, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 45419, 3, 16467,
                                                                       16477, 28082, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 45440, 3, 16477,
                                                                       16487, 28097, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 45461, 3, 16487,
                                                                       16497, 28112, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 45482, 0, 3,
                                                                       45272, 27977, 45293,
                                                                       16517, 16547, 28217,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 45545, 0, 3,
                                                                       45293, 27992, 45314,
                                                                       16547, 16577, 28262,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 45608, 0, 3,
                                                                       45314, 28007, 45335,
                                                                       16577, 16607, 28307,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 45671, 0, 3,
                                                                       45335, 28022, 45356,
                                                                       16607, 16637, 28352,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 45734, 0, 3,
                                                                       45356, 28037, 45377,
                                                                       16637, 16667, 28397,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 45797, 0, 3,
                                                                       45377, 28052, 45398,
                                                                       16667, 16697, 28442,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 45860, 0, 3,
                                                                       45398, 28067, 45419,
                                                                       16697, 16727, 28487,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 45923, 0, 3,
                                                                       45419, 28082, 45440,
                                                                       16727, 16757, 28532,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 45986, 0, 3,
                                                                       45440, 28097, 45461,
                                                                       16757, 16787, 28577,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 46049, 0, 3,
                                                                       45482, 28217, 45545,
                                                                       16847, 16907, 28802,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 46175, 0, 3,
                                                                       45545, 28262, 45608,
                                                                       16907, 16967, 28892,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 46301, 0, 3,
                                                                       45608, 28307, 45671,
                                                                       16967, 17027, 28982,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 46427, 0, 3,
                                                                       45671, 28352, 45734,
                                                                       17027, 17087, 29072,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 46553, 0, 3,
                                                                       45734, 28397, 45797,
                                                                       17087, 17147, 29162,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 46679, 0, 3,
                                                                       45797, 28442, 45860,
                                                                       17147, 17207, 29252,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 46805, 0, 3,
                                                                       45860, 28487, 45923,
                                                                       17207, 17267, 29342,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 46931, 0, 3,
                                                                       45923, 28532, 45986,
                                                                       17267, 17327, 29432,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 47057, 0, 3,
                                                                       46049, 28802, 46175,
                                                                       17447, 17547, 29822,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 47267, 0, 3,
                                                                       46175, 28892, 46301,
                                                                       17547, 17647, 29972,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 47477, 0, 3,
                                                                       46301, 28982, 46427,
                                                                       17647, 17747, 30122,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 47687, 0, 3,
                                                                       46427, 29072, 46553,
                                                                       17747, 17847, 30272,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 47897, 0, 3,
                                                                       46553, 29162, 46679,
                                                                       17847, 17947, 30422,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 48107, 0, 3,
                                                                       46679, 29252, 46805,
                                                                       17947, 18047, 30572,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 48317, 0, 3,
                                                                       46805, 29342, 46931,
                                                                       18047, 18147, 30722,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 48527, 0, 3,
                                                                       47057, 29822, 47267,
                                                                       18347, 18497, 31322,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 48842, 0, 3,
                                                                       47267, 29972, 47477,
                                                                       18497, 18647, 31547,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 49157, 0, 3,
                                                                       47477, 30122, 47687,
                                                                       18647, 18797, 31772,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 49472, 0, 3,
                                                                       47687, 30272, 47897,
                                                                       18797, 18947, 31997,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 49787, 0, 3,
                                                                       47897, 30422, 48107,
                                                                       18947, 19097, 32222,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 50102, 0, 3,
                                                                       48107, 30572, 48317,
                                                                       19097, 19247, 32447,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 50417, 0, 3,
                                                                       48527, 31322, 48842,
                                                                       19547, 19757, 33302,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 50858, 0, 3,
                                                                       48842, 31547, 49157,
                                                                       19757, 19967, 33617,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 51299, 0, 3,
                                                                       49157, 31772, 49472,
                                                                       19967, 20177, 33932,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 51740, 0, 3,
                                                                       49472, 31997, 49787,
                                                                       20177, 20387, 34247,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 52181, 0, 3,
                                                                       49787, 32222, 50102,
                                                                       20387, 20597, 34562,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 52622, 0, 3,
                                                                       50417, 33302, 50858,
                                                                       21017, 21297, 35717,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 53210, 0, 3,
                                                                       50858, 33617, 51299,
                                                                       21297, 21577, 36137,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 53798, 0, 3,
                                                                       51299, 33932, 51740,
                                                                       21577, 21857, 36557,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 54386, 0, 3,
                                                                       51740, 34247, 52181,
                                                                       21857, 22137, 36977,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 54974, 0, 3,
                                                                       52622, 35717, 53210,
                                                                       22697, 23057, 38477,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 55730, 0, 3,
                                                                       53210, 36137, 53798,
                                                                       23057, 23417, 39017,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 56486, 0, 3,
                                                                       53798, 36557, 54386,
                                                                       23417, 23777, 39557,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 57242, 0, 3,
                                                                       54974, 38477, 55730,
                                                                       24497, 24947, 41447,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 58187, 0, 3,
                                                                       55730, 39017, 56486,
                                                                       24947, 25397, 42122,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 59132, 0, 3,
                                                                       57242, 41447, 58187,
                                                                       26297, 26847, 44447,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 60287, 3, 27947,
                                                                       27962, 45272, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 60315, 3, 27962,
                                                                       27977, 45293, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 60343, 3, 27977,
                                                                       27992, 45314, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 60371, 3, 27992,
                                                                       28007, 45335, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 60399, 3, 28007,
                                                                       28022, 45356, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 60427, 3, 28022,
                                                                       28037, 45377, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 60455, 3, 28037,
                                                                       28052, 45398, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 60483, 3, 28052,
                                                                       28067, 45419, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 60511, 3, 28067,
                                                                       28082, 45440, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 60539, 3, 28082,
                                                                       28097, 45461, ncols,
                                                                       gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 60567, 0, 3,
                                                                       60287, 45272, 60315,
                                                                       28127, 28172, 45482,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 60651, 0, 3,
                                                                       60315, 45293, 60343,
                                                                       28172, 28217, 45545,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 60735, 0, 3,
                                                                       60343, 45314, 60371,
                                                                       28217, 28262, 45608,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 60819, 0, 3,
                                                                       60371, 45335, 60399,
                                                                       28262, 28307, 45671,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 60903, 0, 3,
                                                                       60399, 45356, 60427,
                                                                       28307, 28352, 45734,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 60987, 0, 3,
                                                                       60427, 45377, 60455,
                                                                       28352, 28397, 45797,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 61071, 0, 3,
                                                                       60455, 45398, 60483,
                                                                       28397, 28442, 45860,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 61155, 0, 3,
                                                                       60483, 45419, 60511,
                                                                       28442, 28487, 45923,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 61239, 0, 3,
                                                                       60511, 45440, 60539,
                                                                       28487, 28532, 45986,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 61323, 0, 3,
                                                                       60567, 45482, 60651,
                                                                       28622, 28712, 46049,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 61491, 0, 3,
                                                                       60651, 45545, 60735,
                                                                       28712, 28802, 46175,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 61659, 0, 3,
                                                                       60735, 45608, 60819,
                                                                       28802, 28892, 46301,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 61827, 0, 3,
                                                                       60819, 45671, 60903,
                                                                       28892, 28982, 46427,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 61995, 0, 3,
                                                                       60903, 45734, 60987,
                                                                       28982, 29072, 46553,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 62163, 0, 3,
                                                                       60987, 45797, 61071,
                                                                       29072, 29162, 46679,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 62331, 0, 3,
                                                                       61071, 45860, 61155,
                                                                       29162, 29252, 46805,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 62499, 0, 3,
                                                                       61155, 45923, 61239,
                                                                       29252, 29342, 46931,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 62667, 0, 3,
                                                                       61323, 46049, 61491,
                                                                       29522, 29672, 47057,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 62947, 0, 3,
                                                                       61491, 46175, 61659,
                                                                       29672, 29822, 47267,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 63227, 0, 3,
                                                                       61659, 46301, 61827,
                                                                       29822, 29972, 47477,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 63507, 0, 3,
                                                                       61827, 46427, 61995,
                                                                       29972, 30122, 47687,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 63787, 0, 3,
                                                                       61995, 46553, 62163,
                                                                       30122, 30272, 47897,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 64067, 0, 3,
                                                                       62163, 46679, 62331,
                                                                       30272, 30422, 48107,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 64347, 0, 3,
                                                                       62331, 46805, 62499,
                                                                       30422, 30572, 48317,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 64627, 0, 3,
                                                                       62667, 47057, 62947,
                                                                       30872, 31097, 48527,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 65047, 0, 3,
                                                                       62947, 47267, 63227,
                                                                       31097, 31322, 48842,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 65467, 0, 3,
                                                                       63227, 47477, 63507,
                                                                       31322, 31547, 49157,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 65887, 0, 3,
                                                                       63507, 47687, 63787,
                                                                       31547, 31772, 49472,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 66307, 0, 3,
                                                                       63787, 47897, 64067,
                                                                       31772, 31997, 49787,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 66727, 0, 3,
                                                                       64067, 48107, 64347,
                                                                       31997, 32222, 50102,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 67147, 0, 3,
                                                                       64627, 48527, 65047,
                                                                       32672, 32987, 50417,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 67735, 0, 3,
                                                                       65047, 48842, 65467,
                                                                       32987, 33302, 50858,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 68323, 0, 3,
                                                                       65467, 49157, 65887,
                                                                       33302, 33617, 51299,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 68911, 0, 3,
                                                                       65887, 49472, 66307,
                                                                       33617, 33932, 51740,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 69499, 0, 3,
                                                                       66307, 49787, 66727,
                                                                       33932, 34247, 52181,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 70087, 0, 3,
                                                                       67147, 50417, 67735,
                                                                       34877, 35297, 52622,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 70871, 0, 3,
                                                                       67735, 50858, 68323,
                                                                       35297, 35717, 53210,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 71655, 0, 3,
                                                                       68323, 51299, 68911,
                                                                       35717, 36137, 53798,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 72439, 0, 3,
                                                                       68911, 51740, 69499,
                                                                       36137, 36557, 54386,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 73223, 0, 3,
                                                                       70087, 52622, 70871,
                                                                       37397, 37937, 54974,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 74231, 0, 3,
                                                                       70871, 53210, 71655,
                                                                       37937, 38477, 55730,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 75239, 0, 3,
                                                                       71655, 53798, 72439,
                                                                       38477, 39017, 56486,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 76247, 0, 3,
                                                                       73223, 54974, 74231,
                                                                       40097, 40772, 57242,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 77507, 0, 3,
                                                                       74231, 55730, 75239,
                                                                       40772, 41447, 58187,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 78767, 0, 3,
                                                                       76247, 57242, 77507,
                                                                       42797, 43622, 59132,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 80307, 67147, 588, ncols);

                    simdfunc::contract_primitives(buffer, 81168, 70087, 784, ncols);

                    simdfunc::contract_primitives(buffer, 82316, 73223, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 83792, 76247, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 85637, 78767, 1540, ncols);
                }
            }
        }

        simdtrf::transform_i_inner(buffer, 80895, 80307, 21, 1, nmax);

        simdtrf::transform_i_inner(buffer, 81952, 81168, 28, 1, nmax);

        simdtrf::transform_i_inner(buffer, 83324, 82316, 36, 1, nmax);

        simdtrf::transform_i_inner(buffer, 85052, 83792, 45, 1, nmax);

        simdtrf::transform_i_inner(buffer, 87177, 85637, 55, 1, nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 87892, 80895, 81952, 13, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 88711, 81952, 83324, 13, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 89803, 83324, 85052, 13, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 91207, 85052, 87177, 13, nmax);

        simdtrf::compute_hrr_dh(buffer, coordinates, 92962, 87892, 88711, 13, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 94600, 88711, 89803, 13, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 96784, 89803, 91207, 13, nmax);

        simdtrf::compute_hrr_fh(buffer, coordinates, 99592, 92962, 94600, 13, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 102322, 94600, 96784, 13, nmax);

        simdtrf::compute_hrr_gh(buffer, coordinates, 105962, 99592, 102322, 13, nmax);

        simdtrf::transform_h_inner(buffer, 110057, 105962, 15, 13, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 110057, 143, nmax);
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
