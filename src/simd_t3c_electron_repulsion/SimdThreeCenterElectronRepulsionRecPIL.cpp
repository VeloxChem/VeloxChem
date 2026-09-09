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


#include "SimdThreeCenterElectronRepulsionRecPIL.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSKD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKS.hpp"
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
#include "SimdTransferPI.hpp"
#include "SimdTransformI.hpp"
#include "SimdTransformL.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_pil_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_pil_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 91716, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 663 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 91716, 85657, 3356, dimensions);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1297, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1300, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1303, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1306, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1309, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1312, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1315, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1318, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1321, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1324, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1327, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1330, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1333, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1336, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1339, 3, 9, 29,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1348, 3, 10, 32,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1357, 3, 11, 35,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1366, 3, 12, 38,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1375, 3, 13, 41,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1384, 3, 14, 44,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1393, 3, 15, 47,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1402, 3, 16, 50,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1411, 3, 17, 53,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1420, 3, 18, 56,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1429, 3, 19, 59,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1438, 3, 20, 62,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1447, 3, 21, 65,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1456, 3, 29, 80,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1474, 3, 32, 86,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1492, 3, 35, 92,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1510, 3, 38, 98,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1528, 3, 41, 104,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1546, 3, 44, 110,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1564, 3, 47, 116,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1582, 3, 50, 122,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1600, 3, 53, 128,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1618, 3, 56, 134,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1636, 3, 59, 140,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1654, 3, 62, 146,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1672, 3, 80, 172,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1702, 3, 86, 182,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1732, 3, 92, 192,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1762, 3, 98, 202,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1792, 3, 104, 212,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1822, 3, 110, 222,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1852, 3, 116, 232,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1882, 3, 122, 242,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1912, 3, 128, 252,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1942, 3, 134, 262,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1972, 3, 140, 272,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2002, 3, 172, 312,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2047, 3, 182, 327,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2092, 3, 192, 342,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2137, 3, 202, 357,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2182, 3, 212, 372,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2227, 3, 222, 387,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2272, 3, 232, 402,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2317, 3, 242, 417,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2362, 3, 252, 432,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2407, 3, 262, 447,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2452, 3, 312, 504,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2515, 3, 327, 525,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2578, 3, 342, 546,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2641, 3, 357, 567,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2704, 3, 372, 588,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2767, 3, 387, 609,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2830, 3, 402, 630,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2893, 3, 417, 651,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2956, 3, 432, 672,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3019, 3, 504, 749,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3103, 3, 525, 777,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3187, 3, 546, 805,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3271, 3, 567, 833,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3355, 3, 588, 861,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3439, 3, 609, 889,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3523, 3, 630, 917,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3607, 3, 651, 945,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 3691, 3, 749,
                                                                       1045, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 3799, 3, 777,
                                                                       1081, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 3907, 3, 805,
                                                                       1117, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4015, 3, 833,
                                                                       1153, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4123, 3, 861,
                                                                       1189, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4231, 3, 889,
                                                                       1225, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4339, 3, 917,
                                                                       1261, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4447, 3, 7, 8,
                                                                       1297, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4453, 3, 8, 9,
                                                                       1300, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4459, 3, 9, 10,
                                                                       1303, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4465, 3, 10, 11,
                                                                       1306, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4471, 3, 11, 12,
                                                                       1309, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4477, 3, 12, 13,
                                                                       1312, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4483, 3, 13, 14,
                                                                       1315, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4489, 3, 14, 15,
                                                                       1318, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4495, 3, 15, 16,
                                                                       1321, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4501, 3, 16, 17,
                                                                       1324, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4507, 3, 17, 18,
                                                                       1327, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4513, 3, 18, 19,
                                                                       1330, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4519, 3, 19, 20,
                                                                       1333, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4525, 3, 20, 21,
                                                                       1336, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4531, 0, 3, 4447,
                                                                       1297, 4453, 1339, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4549, 0, 3, 4453,
                                                                       1300, 4459, 1348, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4567, 0, 3, 4459,
                                                                       1303, 4465, 1357, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4585, 0, 3, 4465,
                                                                       1306, 4471, 1366, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4603, 0, 3, 4471,
                                                                       1309, 4477, 1375, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4621, 0, 3, 4477,
                                                                       1312, 4483, 1384, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4639, 0, 3, 4483,
                                                                       1315, 4489, 1393, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4657, 0, 3, 4489,
                                                                       1318, 4495, 1402, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4675, 0, 3, 4495,
                                                                       1321, 4501, 1411, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4693, 0, 3, 4501,
                                                                       1324, 4507, 1420, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4711, 0, 3, 4507,
                                                                       1327, 4513, 1429, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4729, 0, 3, 4513,
                                                                       1330, 4519, 1438, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4747, 0, 3, 4519,
                                                                       1333, 4525, 1447, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4765, 0, 3, 4531,
                                                                       1339, 4549, 68, 74, 1456,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4801, 0, 3, 4549,
                                                                       1348, 4567, 74, 80, 1474,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4837, 0, 3, 4567,
                                                                       1357, 4585, 80, 86, 1492,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4873, 0, 3, 4585,
                                                                       1366, 4603, 86, 92, 1510,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4909, 0, 3, 4603,
                                                                       1375, 4621, 92, 98, 1528,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4945, 0, 3, 4621,
                                                                       1384, 4639, 98, 104, 1546,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4981, 0, 3, 4639,
                                                                       1393, 4657, 104, 110,
                                                                       1564, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5017, 0, 3, 4657,
                                                                       1402, 4675, 110, 116,
                                                                       1582, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5053, 0, 3, 4675,
                                                                       1411, 4693, 116, 122,
                                                                       1600, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5089, 0, 3, 4693,
                                                                       1420, 4711, 122, 128,
                                                                       1618, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5125, 0, 3, 4711,
                                                                       1429, 4729, 128, 134,
                                                                       1636, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5161, 0, 3, 4729,
                                                                       1438, 4747, 134, 140,
                                                                       1654, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5197, 0, 3, 4765,
                                                                       1456, 4801, 152, 162,
                                                                       1672, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5257, 0, 3, 4801,
                                                                       1474, 4837, 162, 172,
                                                                       1702, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5317, 0, 3, 4837,
                                                                       1492, 4873, 172, 182,
                                                                       1732, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5377, 0, 3, 4873,
                                                                       1510, 4909, 182, 192,
                                                                       1762, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5437, 0, 3, 4909,
                                                                       1528, 4945, 192, 202,
                                                                       1792, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5497, 0, 3, 4945,
                                                                       1546, 4981, 202, 212,
                                                                       1822, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5557, 0, 3, 4981,
                                                                       1564, 5017, 212, 222,
                                                                       1852, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5617, 0, 3, 5017,
                                                                       1582, 5053, 222, 232,
                                                                       1882, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5677, 0, 3, 5053,
                                                                       1600, 5089, 232, 242,
                                                                       1912, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5737, 0, 3, 5089,
                                                                       1618, 5125, 242, 252,
                                                                       1942, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5797, 0, 3, 5125,
                                                                       1636, 5161, 252, 262,
                                                                       1972, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5857, 0, 3, 5197,
                                                                       1672, 5257, 282, 297,
                                                                       2002, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5947, 0, 3, 5257,
                                                                       1702, 5317, 297, 312,
                                                                       2047, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6037, 0, 3, 5317,
                                                                       1732, 5377, 312, 327,
                                                                       2092, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6127, 0, 3, 5377,
                                                                       1762, 5437, 327, 342,
                                                                       2137, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6217, 0, 3, 5437,
                                                                       1792, 5497, 342, 357,
                                                                       2182, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6307, 0, 3, 5497,
                                                                       1822, 5557, 357, 372,
                                                                       2227, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6397, 0, 3, 5557,
                                                                       1852, 5617, 372, 387,
                                                                       2272, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6487, 0, 3, 5617,
                                                                       1882, 5677, 387, 402,
                                                                       2317, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6577, 0, 3, 5677,
                                                                       1912, 5737, 402, 417,
                                                                       2362, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6667, 0, 3, 5737,
                                                                       1942, 5797, 417, 432,
                                                                       2407, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 6757, 0, 3, 5857,
                                                                       2002, 5947, 462, 483,
                                                                       2452, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 6883, 0, 3, 5947,
                                                                       2047, 6037, 483, 504,
                                                                       2515, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 7009, 0, 3, 6037,
                                                                       2092, 6127, 504, 525,
                                                                       2578, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 7135, 0, 3, 6127,
                                                                       2137, 6217, 525, 546,
                                                                       2641, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 7261, 0, 3, 6217,
                                                                       2182, 6307, 546, 567,
                                                                       2704, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 7387, 0, 3, 6307,
                                                                       2227, 6397, 567, 588,
                                                                       2767, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 7513, 0, 3, 6397,
                                                                       2272, 6487, 588, 609,
                                                                       2830, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 7639, 0, 3, 6487,
                                                                       2317, 6577, 609, 630,
                                                                       2893, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 7765, 0, 3, 6577,
                                                                       2362, 6667, 630, 651,
                                                                       2956, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 7891, 0, 3, 6757,
                                                                       2452, 6883, 693, 721,
                                                                       3019, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 8059, 0, 3, 6883,
                                                                       2515, 7009, 721, 749,
                                                                       3103, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 8227, 0, 3, 7009,
                                                                       2578, 7135, 749, 777,
                                                                       3187, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 8395, 0, 3, 7135,
                                                                       2641, 7261, 777, 805,
                                                                       3271, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 8563, 0, 3, 7261,
                                                                       2704, 7387, 805, 833,
                                                                       3355, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 8731, 0, 3, 7387,
                                                                       2767, 7513, 833, 861,
                                                                       3439, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 8899, 0, 3, 7513,
                                                                       2830, 7639, 861, 889,
                                                                       3523, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 9067, 0, 3, 7639,
                                                                       2893, 7765, 889, 917,
                                                                       3607, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 9235, 0, 3, 7891,
                                                                       3019, 8059, 973, 1009,
                                                                       3691, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 9451, 0, 3, 8059,
                                                                       3103, 8227, 1009, 1045,
                                                                       3799, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 9667, 0, 3, 8227,
                                                                       3187, 8395, 1045, 1081,
                                                                       3907, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 9883, 0, 3, 8395,
                                                                       3271, 8563, 1081, 1117,
                                                                       4015, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 10099, 0, 3, 8563,
                                                                       3355, 8731, 1117, 1153,
                                                                       4123, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 10315, 0, 3, 8731,
                                                                       3439, 8899, 1153, 1189,
                                                                       4231, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 10531, 0, 3, 8899,
                                                                       3523, 9067, 1189, 1225,
                                                                       4339, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 10747, 3, 1297,
                                                                       1300, 4459, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 10757, 3, 1300,
                                                                       1303, 4465, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 10767, 3, 1303,
                                                                       1306, 4471, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 10777, 3, 1306,
                                                                       1309, 4477, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 10787, 3, 1309,
                                                                       1312, 4483, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 10797, 3, 1312,
                                                                       1315, 4489, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 10807, 3, 1315,
                                                                       1318, 4495, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 10817, 3, 1318,
                                                                       1321, 4501, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 10827, 3, 1321,
                                                                       1324, 4507, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 10837, 3, 1324,
                                                                       1327, 4513, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 10847, 3, 1327,
                                                                       1330, 4519, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 10857, 3, 1330,
                                                                       1333, 4525, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 10867, 0, 3,
                                                                       10747, 4459, 10757, 4567,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 10897, 0, 3,
                                                                       10757, 4465, 10767, 4585,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 10927, 0, 3,
                                                                       10767, 4471, 10777, 4603,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 10957, 0, 3,
                                                                       10777, 4477, 10787, 4621,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 10987, 0, 3,
                                                                       10787, 4483, 10797, 4639,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 11017, 0, 3,
                                                                       10797, 4489, 10807, 4657,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 11047, 0, 3,
                                                                       10807, 4495, 10817, 4675,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 11077, 0, 3,
                                                                       10817, 4501, 10827, 4693,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 11107, 0, 3,
                                                                       10827, 4507, 10837, 4711,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 11137, 0, 3,
                                                                       10837, 4513, 10847, 4729,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 11167, 0, 3,
                                                                       10847, 4519, 10857, 4747,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 11197, 0, 3,
                                                                       10867, 4567, 10897, 1456,
                                                                       1474, 4837, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 11257, 0, 3,
                                                                       10897, 4585, 10927, 1474,
                                                                       1492, 4873, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 11317, 0, 3,
                                                                       10927, 4603, 10957, 1492,
                                                                       1510, 4909, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 11377, 0, 3,
                                                                       10957, 4621, 10987, 1510,
                                                                       1528, 4945, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 11437, 0, 3,
                                                                       10987, 4639, 11017, 1528,
                                                                       1546, 4981, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 11497, 0, 3,
                                                                       11017, 4657, 11047, 1546,
                                                                       1564, 5017, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 11557, 0, 3,
                                                                       11047, 4675, 11077, 1564,
                                                                       1582, 5053, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 11617, 0, 3,
                                                                       11077, 4693, 11107, 1582,
                                                                       1600, 5089, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 11677, 0, 3,
                                                                       11107, 4711, 11137, 1600,
                                                                       1618, 5125, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 11737, 0, 3,
                                                                       11137, 4729, 11167, 1618,
                                                                       1636, 5161, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 11797, 0, 3,
                                                                       11197, 4837, 11257, 1672,
                                                                       1702, 5317, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 11897, 0, 3,
                                                                       11257, 4873, 11317, 1702,
                                                                       1732, 5377, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 11997, 0, 3,
                                                                       11317, 4909, 11377, 1732,
                                                                       1762, 5437, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 12097, 0, 3,
                                                                       11377, 4945, 11437, 1762,
                                                                       1792, 5497, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 12197, 0, 3,
                                                                       11437, 4981, 11497, 1792,
                                                                       1822, 5557, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 12297, 0, 3,
                                                                       11497, 5017, 11557, 1822,
                                                                       1852, 5617, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 12397, 0, 3,
                                                                       11557, 5053, 11617, 1852,
                                                                       1882, 5677, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 12497, 0, 3,
                                                                       11617, 5089, 11677, 1882,
                                                                       1912, 5737, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 12597, 0, 3,
                                                                       11677, 5125, 11737, 1912,
                                                                       1942, 5797, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 12697, 0, 3,
                                                                       11797, 5317, 11897, 2002,
                                                                       2047, 6037, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 12847, 0, 3,
                                                                       11897, 5377, 11997, 2047,
                                                                       2092, 6127, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 12997, 0, 3,
                                                                       11997, 5437, 12097, 2092,
                                                                       2137, 6217, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 13147, 0, 3,
                                                                       12097, 5497, 12197, 2137,
                                                                       2182, 6307, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 13297, 0, 3,
                                                                       12197, 5557, 12297, 2182,
                                                                       2227, 6397, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 13447, 0, 3,
                                                                       12297, 5617, 12397, 2227,
                                                                       2272, 6487, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 13597, 0, 3,
                                                                       12397, 5677, 12497, 2272,
                                                                       2317, 6577, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 13747, 0, 3,
                                                                       12497, 5737, 12597, 2317,
                                                                       2362, 6667, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 13897, 0, 3,
                                                                       12697, 6037, 12847, 2452,
                                                                       2515, 7009, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 14107, 0, 3,
                                                                       12847, 6127, 12997, 2515,
                                                                       2578, 7135, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 14317, 0, 3,
                                                                       12997, 6217, 13147, 2578,
                                                                       2641, 7261, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 14527, 0, 3,
                                                                       13147, 6307, 13297, 2641,
                                                                       2704, 7387, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 14737, 0, 3,
                                                                       13297, 6397, 13447, 2704,
                                                                       2767, 7513, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 14947, 0, 3,
                                                                       13447, 6487, 13597, 2767,
                                                                       2830, 7639, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 15157, 0, 3,
                                                                       13597, 6577, 13747, 2830,
                                                                       2893, 7765, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 15367, 0, 3,
                                                                       13897, 7009, 14107, 3019,
                                                                       3103, 8227, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 15647, 0, 3,
                                                                       14107, 7135, 14317, 3103,
                                                                       3187, 8395, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 15927, 0, 3,
                                                                       14317, 7261, 14527, 3187,
                                                                       3271, 8563, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 16207, 0, 3,
                                                                       14527, 7387, 14737, 3271,
                                                                       3355, 8731, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 16487, 0, 3,
                                                                       14737, 7513, 14947, 3355,
                                                                       3439, 8899, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 16767, 0, 3,
                                                                       14947, 7639, 15157, 3439,
                                                                       3523, 9067, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 17047, 0, 3,
                                                                       15367, 8227, 15647, 3691,
                                                                       3799, 9667, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 17407, 0, 3,
                                                                       15647, 8395, 15927, 3799,
                                                                       3907, 9883, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 17767, 0, 3,
                                                                       15927, 8563, 16207, 3907,
                                                                       4015, 10099, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 18127, 0, 3,
                                                                       16207, 8731, 16487, 4015,
                                                                       4123, 10315, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 18487, 0, 3,
                                                                       16487, 8899, 16767, 4123,
                                                                       4231, 10531, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18847, 3, 4447,
                                                                       4453, 10747, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18862, 3, 4453,
                                                                       4459, 10757, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18877, 3, 4459,
                                                                       4465, 10767, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18892, 3, 4465,
                                                                       4471, 10777, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18907, 3, 4471,
                                                                       4477, 10787, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18922, 3, 4477,
                                                                       4483, 10797, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18937, 3, 4483,
                                                                       4489, 10807, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18952, 3, 4489,
                                                                       4495, 10817, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18967, 3, 4495,
                                                                       4501, 10827, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18982, 3, 4501,
                                                                       4507, 10837, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18997, 3, 4507,
                                                                       4513, 10847, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 19012, 3, 4513,
                                                                       4519, 10857, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19027, 0, 3,
                                                                       18847, 10747, 18862, 4531,
                                                                       4549, 10867, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19072, 0, 3,
                                                                       18862, 10757, 18877, 4549,
                                                                       4567, 10897, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19117, 0, 3,
                                                                       18877, 10767, 18892, 4567,
                                                                       4585, 10927, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19162, 0, 3,
                                                                       18892, 10777, 18907, 4585,
                                                                       4603, 10957, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19207, 0, 3,
                                                                       18907, 10787, 18922, 4603,
                                                                       4621, 10987, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19252, 0, 3,
                                                                       18922, 10797, 18937, 4621,
                                                                       4639, 11017, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19297, 0, 3,
                                                                       18937, 10807, 18952, 4639,
                                                                       4657, 11047, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19342, 0, 3,
                                                                       18952, 10817, 18967, 4657,
                                                                       4675, 11077, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19387, 0, 3,
                                                                       18967, 10827, 18982, 4675,
                                                                       4693, 11107, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19432, 0, 3,
                                                                       18982, 10837, 18997, 4693,
                                                                       4711, 11137, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19477, 0, 3,
                                                                       18997, 10847, 19012, 4711,
                                                                       4729, 11167, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 19522, 0, 3,
                                                                       19027, 10867, 19072, 4765,
                                                                       4801, 11197, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 19612, 0, 3,
                                                                       19072, 10897, 19117, 4801,
                                                                       4837, 11257, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 19702, 0, 3,
                                                                       19117, 10927, 19162, 4837,
                                                                       4873, 11317, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 19792, 0, 3,
                                                                       19162, 10957, 19207, 4873,
                                                                       4909, 11377, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 19882, 0, 3,
                                                                       19207, 10987, 19252, 4909,
                                                                       4945, 11437, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 19972, 0, 3,
                                                                       19252, 11017, 19297, 4945,
                                                                       4981, 11497, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 20062, 0, 3,
                                                                       19297, 11047, 19342, 4981,
                                                                       5017, 11557, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 20152, 0, 3,
                                                                       19342, 11077, 19387, 5017,
                                                                       5053, 11617, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 20242, 0, 3,
                                                                       19387, 11107, 19432, 5053,
                                                                       5089, 11677, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 20332, 0, 3,
                                                                       19432, 11137, 19477, 5089,
                                                                       5125, 11737, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 20422, 0, 3,
                                                                       19522, 11197, 19612, 5197,
                                                                       5257, 11797, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 20572, 0, 3,
                                                                       19612, 11257, 19702, 5257,
                                                                       5317, 11897, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 20722, 0, 3,
                                                                       19702, 11317, 19792, 5317,
                                                                       5377, 11997, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 20872, 0, 3,
                                                                       19792, 11377, 19882, 5377,
                                                                       5437, 12097, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 21022, 0, 3,
                                                                       19882, 11437, 19972, 5437,
                                                                       5497, 12197, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 21172, 0, 3,
                                                                       19972, 11497, 20062, 5497,
                                                                       5557, 12297, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 21322, 0, 3,
                                                                       20062, 11557, 20152, 5557,
                                                                       5617, 12397, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 21472, 0, 3,
                                                                       20152, 11617, 20242, 5617,
                                                                       5677, 12497, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 21622, 0, 3,
                                                                       20242, 11677, 20332, 5677,
                                                                       5737, 12597, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 21772, 0, 3,
                                                                       20422, 11797, 20572, 5857,
                                                                       5947, 12697, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 21997, 0, 3,
                                                                       20572, 11897, 20722, 5947,
                                                                       6037, 12847, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 22222, 0, 3,
                                                                       20722, 11997, 20872, 6037,
                                                                       6127, 12997, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 22447, 0, 3,
                                                                       20872, 12097, 21022, 6127,
                                                                       6217, 13147, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 22672, 0, 3,
                                                                       21022, 12197, 21172, 6217,
                                                                       6307, 13297, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 22897, 0, 3,
                                                                       21172, 12297, 21322, 6307,
                                                                       6397, 13447, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 23122, 0, 3,
                                                                       21322, 12397, 21472, 6397,
                                                                       6487, 13597, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 23347, 0, 3,
                                                                       21472, 12497, 21622, 6487,
                                                                       6577, 13747, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 23572, 0, 3,
                                                                       21772, 12697, 21997, 6757,
                                                                       6883, 13897, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 23887, 0, 3,
                                                                       21997, 12847, 22222, 6883,
                                                                       7009, 14107, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 24202, 0, 3,
                                                                       22222, 12997, 22447, 7009,
                                                                       7135, 14317, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 24517, 0, 3,
                                                                       22447, 13147, 22672, 7135,
                                                                       7261, 14527, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 24832, 0, 3,
                                                                       22672, 13297, 22897, 7261,
                                                                       7387, 14737, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 25147, 0, 3,
                                                                       22897, 13447, 23122, 7387,
                                                                       7513, 14947, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 25462, 0, 3,
                                                                       23122, 13597, 23347, 7513,
                                                                       7639, 15157, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 25777, 0, 3,
                                                                       23572, 13897, 23887, 7891,
                                                                       8059, 15367, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 26197, 0, 3,
                                                                       23887, 14107, 24202, 8059,
                                                                       8227, 15647, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 26617, 0, 3,
                                                                       24202, 14317, 24517, 8227,
                                                                       8395, 15927, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 27037, 0, 3,
                                                                       24517, 14527, 24832, 8395,
                                                                       8563, 16207, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 27457, 0, 3,
                                                                       24832, 14737, 25147, 8563,
                                                                       8731, 16487, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 27877, 0, 3,
                                                                       25147, 14947, 25462, 8731,
                                                                       8899, 16767, ncols, gamma,
                                                                       p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 28297, 0, 3,
                                                                       25777, 15367, 26197, 9235,
                                                                       9451, 17047, ncols, gamma,
                                                                       p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 28837, 0, 3,
                                                                       26197, 15647, 26617, 9451,
                                                                       9667, 17407, ncols, gamma,
                                                                       p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 29377, 0, 3,
                                                                       26617, 15927, 27037, 9667,
                                                                       9883, 17767, ncols, gamma,
                                                                       p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 29917, 0, 3,
                                                                       27037, 16207, 27457, 9883,
                                                                       10099, 18127, ncols,
                                                                       gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 30457, 0, 3,
                                                                       27457, 16487, 27877,
                                                                       10099, 10315, 18487,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 30997, 3, 10747,
                                                                       10757, 18877, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 31018, 3, 10757,
                                                                       10767, 18892, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 31039, 3, 10767,
                                                                       10777, 18907, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 31060, 3, 10777,
                                                                       10787, 18922, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 31081, 3, 10787,
                                                                       10797, 18937, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 31102, 3, 10797,
                                                                       10807, 18952, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 31123, 3, 10807,
                                                                       10817, 18967, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 31144, 3, 10817,
                                                                       10827, 18982, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 31165, 3, 10827,
                                                                       10837, 18997, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 31186, 3, 10837,
                                                                       10847, 19012, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 31207, 0, 3,
                                                                       30997, 18877, 31018,
                                                                       10867, 10897, 19117,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 31270, 0, 3,
                                                                       31018, 18892, 31039,
                                                                       10897, 10927, 19162,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 31333, 0, 3,
                                                                       31039, 18907, 31060,
                                                                       10927, 10957, 19207,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 31396, 0, 3,
                                                                       31060, 18922, 31081,
                                                                       10957, 10987, 19252,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 31459, 0, 3,
                                                                       31081, 18937, 31102,
                                                                       10987, 11017, 19297,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 31522, 0, 3,
                                                                       31102, 18952, 31123,
                                                                       11017, 11047, 19342,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 31585, 0, 3,
                                                                       31123, 18967, 31144,
                                                                       11047, 11077, 19387,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 31648, 0, 3,
                                                                       31144, 18982, 31165,
                                                                       11077, 11107, 19432,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 31711, 0, 3,
                                                                       31165, 18997, 31186,
                                                                       11107, 11137, 19477,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 31774, 0, 3,
                                                                       31207, 19117, 31270,
                                                                       11197, 11257, 19702,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 31900, 0, 3,
                                                                       31270, 19162, 31333,
                                                                       11257, 11317, 19792,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 32026, 0, 3,
                                                                       31333, 19207, 31396,
                                                                       11317, 11377, 19882,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 32152, 0, 3,
                                                                       31396, 19252, 31459,
                                                                       11377, 11437, 19972,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 32278, 0, 3,
                                                                       31459, 19297, 31522,
                                                                       11437, 11497, 20062,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 32404, 0, 3,
                                                                       31522, 19342, 31585,
                                                                       11497, 11557, 20152,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 32530, 0, 3,
                                                                       31585, 19387, 31648,
                                                                       11557, 11617, 20242,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 32656, 0, 3,
                                                                       31648, 19432, 31711,
                                                                       11617, 11677, 20332,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 32782, 0, 3,
                                                                       31774, 19702, 31900,
                                                                       11797, 11897, 20722,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 32992, 0, 3,
                                                                       31900, 19792, 32026,
                                                                       11897, 11997, 20872,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 33202, 0, 3,
                                                                       32026, 19882, 32152,
                                                                       11997, 12097, 21022,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 33412, 0, 3,
                                                                       32152, 19972, 32278,
                                                                       12097, 12197, 21172,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 33622, 0, 3,
                                                                       32278, 20062, 32404,
                                                                       12197, 12297, 21322,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 33832, 0, 3,
                                                                       32404, 20152, 32530,
                                                                       12297, 12397, 21472,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 34042, 0, 3,
                                                                       32530, 20242, 32656,
                                                                       12397, 12497, 21622,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 34252, 0, 3,
                                                                       32782, 20722, 32992,
                                                                       12697, 12847, 22222,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 34567, 0, 3,
                                                                       32992, 20872, 33202,
                                                                       12847, 12997, 22447,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 34882, 0, 3,
                                                                       33202, 21022, 33412,
                                                                       12997, 13147, 22672,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 35197, 0, 3,
                                                                       33412, 21172, 33622,
                                                                       13147, 13297, 22897,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 35512, 0, 3,
                                                                       33622, 21322, 33832,
                                                                       13297, 13447, 23122,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 35827, 0, 3,
                                                                       33832, 21472, 34042,
                                                                       13447, 13597, 23347,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 36142, 0, 3,
                                                                       34252, 22222, 34567,
                                                                       13897, 14107, 24202,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 36583, 0, 3,
                                                                       34567, 22447, 34882,
                                                                       14107, 14317, 24517,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 37024, 0, 3,
                                                                       34882, 22672, 35197,
                                                                       14317, 14527, 24832,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 37465, 0, 3,
                                                                       35197, 22897, 35512,
                                                                       14527, 14737, 25147,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 37906, 0, 3,
                                                                       35512, 23122, 35827,
                                                                       14737, 14947, 25462,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 38347, 0, 3,
                                                                       36142, 24202, 36583,
                                                                       15367, 15647, 26617,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 38935, 0, 3,
                                                                       36583, 24517, 37024,
                                                                       15647, 15927, 27037,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 39523, 0, 3,
                                                                       37024, 24832, 37465,
                                                                       15927, 16207, 27457,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 40111, 0, 3,
                                                                       37465, 25147, 37906,
                                                                       16207, 16487, 27877,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 40699, 0, 3,
                                                                       38347, 26617, 38935,
                                                                       17047, 17407, 29377,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 41455, 0, 3,
                                                                       38935, 27037, 39523,
                                                                       17407, 17767, 29917,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 42211, 0, 3,
                                                                       39523, 27457, 40111,
                                                                       17767, 18127, 30457,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 42967, 3, 18847,
                                                                       18862, 30997, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 42995, 3, 18862,
                                                                       18877, 31018, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 43023, 3, 18877,
                                                                       18892, 31039, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 43051, 3, 18892,
                                                                       18907, 31060, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 43079, 3, 18907,
                                                                       18922, 31081, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 43107, 3, 18922,
                                                                       18937, 31102, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 43135, 3, 18937,
                                                                       18952, 31123, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 43163, 3, 18952,
                                                                       18967, 31144, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 43191, 3, 18967,
                                                                       18982, 31165, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 43219, 3, 18982,
                                                                       18997, 31186, ncols,
                                                                       gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 43247, 0, 3,
                                                                       42967, 30997, 42995,
                                                                       19027, 19072, 31207,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 43331, 0, 3,
                                                                       42995, 31018, 43023,
                                                                       19072, 19117, 31270,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 43415, 0, 3,
                                                                       43023, 31039, 43051,
                                                                       19117, 19162, 31333,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 43499, 0, 3,
                                                                       43051, 31060, 43079,
                                                                       19162, 19207, 31396,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 43583, 0, 3,
                                                                       43079, 31081, 43107,
                                                                       19207, 19252, 31459,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 43667, 0, 3,
                                                                       43107, 31102, 43135,
                                                                       19252, 19297, 31522,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 43751, 0, 3,
                                                                       43135, 31123, 43163,
                                                                       19297, 19342, 31585,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 43835, 0, 3,
                                                                       43163, 31144, 43191,
                                                                       19342, 19387, 31648,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 43919, 0, 3,
                                                                       43191, 31165, 43219,
                                                                       19387, 19432, 31711,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 44003, 0, 3,
                                                                       43247, 31207, 43331,
                                                                       19522, 19612, 31774,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 44171, 0, 3,
                                                                       43331, 31270, 43415,
                                                                       19612, 19702, 31900,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 44339, 0, 3,
                                                                       43415, 31333, 43499,
                                                                       19702, 19792, 32026,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 44507, 0, 3,
                                                                       43499, 31396, 43583,
                                                                       19792, 19882, 32152,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 44675, 0, 3,
                                                                       43583, 31459, 43667,
                                                                       19882, 19972, 32278,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 44843, 0, 3,
                                                                       43667, 31522, 43751,
                                                                       19972, 20062, 32404,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 45011, 0, 3,
                                                                       43751, 31585, 43835,
                                                                       20062, 20152, 32530,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 45179, 0, 3,
                                                                       43835, 31648, 43919,
                                                                       20152, 20242, 32656,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 45347, 0, 3,
                                                                       44003, 31774, 44171,
                                                                       20422, 20572, 32782,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 45627, 0, 3,
                                                                       44171, 31900, 44339,
                                                                       20572, 20722, 32992,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 45907, 0, 3,
                                                                       44339, 32026, 44507,
                                                                       20722, 20872, 33202,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 46187, 0, 3,
                                                                       44507, 32152, 44675,
                                                                       20872, 21022, 33412,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 46467, 0, 3,
                                                                       44675, 32278, 44843,
                                                                       21022, 21172, 33622,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 46747, 0, 3,
                                                                       44843, 32404, 45011,
                                                                       21172, 21322, 33832,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 47027, 0, 3,
                                                                       45011, 32530, 45179,
                                                                       21322, 21472, 34042,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 47307, 0, 3,
                                                                       45347, 32782, 45627,
                                                                       21772, 21997, 34252,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 47727, 0, 3,
                                                                       45627, 32992, 45907,
                                                                       21997, 22222, 34567,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 48147, 0, 3,
                                                                       45907, 33202, 46187,
                                                                       22222, 22447, 34882,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 48567, 0, 3,
                                                                       46187, 33412, 46467,
                                                                       22447, 22672, 35197,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 48987, 0, 3,
                                                                       46467, 33622, 46747,
                                                                       22672, 22897, 35512,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 49407, 0, 3,
                                                                       46747, 33832, 47027,
                                                                       22897, 23122, 35827,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 49827, 0, 3,
                                                                       47307, 34252, 47727,
                                                                       23572, 23887, 36142,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 50415, 0, 3,
                                                                       47727, 34567, 48147,
                                                                       23887, 24202, 36583,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 51003, 0, 3,
                                                                       48147, 34882, 48567,
                                                                       24202, 24517, 37024,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 51591, 0, 3,
                                                                       48567, 35197, 48987,
                                                                       24517, 24832, 37465,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 52179, 0, 3,
                                                                       48987, 35512, 49407,
                                                                       24832, 25147, 37906,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 52767, 0, 3,
                                                                       49827, 36142, 50415,
                                                                       25777, 26197, 38347,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 53551, 0, 3,
                                                                       50415, 36583, 51003,
                                                                       26197, 26617, 38935,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 54335, 0, 3,
                                                                       51003, 37024, 51591,
                                                                       26617, 27037, 39523,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 55119, 0, 3,
                                                                       51591, 37465, 52179,
                                                                       27037, 27457, 40111,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 55903, 0, 3,
                                                                       52767, 38347, 53551,
                                                                       28297, 28837, 40699,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 56911, 0, 3,
                                                                       53551, 38935, 54335,
                                                                       28837, 29377, 41455,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 57919, 0, 3,
                                                                       54335, 39523, 55119,
                                                                       29377, 29917, 42211,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 58927, 3, 30997,
                                                                       31018, 43023, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 58963, 3, 31018,
                                                                       31039, 43051, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 58999, 3, 31039,
                                                                       31060, 43079, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 59035, 3, 31060,
                                                                       31081, 43107, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 59071, 3, 31081,
                                                                       31102, 43135, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 59107, 3, 31102,
                                                                       31123, 43163, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 59143, 3, 31123,
                                                                       31144, 43191, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 59179, 3, 31144,
                                                                       31165, 43219, ncols,
                                                                       gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 59215, 0, 3,
                                                                       58927, 43023, 58963,
                                                                       31207, 31270, 43415,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 59323, 0, 3,
                                                                       58963, 43051, 58999,
                                                                       31270, 31333, 43499,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 59431, 0, 3,
                                                                       58999, 43079, 59035,
                                                                       31333, 31396, 43583,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 59539, 0, 3,
                                                                       59035, 43107, 59071,
                                                                       31396, 31459, 43667,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 59647, 0, 3,
                                                                       59071, 43135, 59107,
                                                                       31459, 31522, 43751,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 59755, 0, 3,
                                                                       59107, 43163, 59143,
                                                                       31522, 31585, 43835,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 59863, 0, 3,
                                                                       59143, 43191, 59179,
                                                                       31585, 31648, 43919,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 59971, 0, 3,
                                                                       59215, 43415, 59323,
                                                                       31774, 31900, 44339,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 60187, 0, 3,
                                                                       59323, 43499, 59431,
                                                                       31900, 32026, 44507,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 60403, 0, 3,
                                                                       59431, 43583, 59539,
                                                                       32026, 32152, 44675,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 60619, 0, 3,
                                                                       59539, 43667, 59647,
                                                                       32152, 32278, 44843,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 60835, 0, 3,
                                                                       59647, 43751, 59755,
                                                                       32278, 32404, 45011,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 61051, 0, 3,
                                                                       59755, 43835, 59863,
                                                                       32404, 32530, 45179,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 61267, 0, 3,
                                                                       59971, 44339, 60187,
                                                                       32782, 32992, 45907,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 61627, 0, 3,
                                                                       60187, 44507, 60403,
                                                                       32992, 33202, 46187,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 61987, 0, 3,
                                                                       60403, 44675, 60619,
                                                                       33202, 33412, 46467,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 62347, 0, 3,
                                                                       60619, 44843, 60835,
                                                                       33412, 33622, 46747,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 62707, 0, 3,
                                                                       60835, 45011, 61051,
                                                                       33622, 33832, 47027,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 63067, 0, 3,
                                                                       61267, 45907, 61627,
                                                                       34252, 34567, 48147,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 63607, 0, 3,
                                                                       61627, 46187, 61987,
                                                                       34567, 34882, 48567,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 64147, 0, 3,
                                                                       61987, 46467, 62347,
                                                                       34882, 35197, 48987,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 64687, 0, 3,
                                                                       62347, 46747, 62707,
                                                                       35197, 35512, 49407,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 65227, 0, 3,
                                                                       63067, 48147, 63607,
                                                                       36142, 36583, 51003,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 65983, 0, 3,
                                                                       63607, 48567, 64147,
                                                                       36583, 37024, 51591,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 66739, 0, 3,
                                                                       64147, 48987, 64687,
                                                                       37024, 37465, 52179,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 67495, 0, 3,
                                                                       65227, 51003, 65983,
                                                                       38347, 38935, 54335,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 68503, 0, 3,
                                                                       65983, 51591, 66739,
                                                                       38935, 39523, 55119,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 69511, 0, 3,
                                                                       67495, 54335, 68503,
                                                                       40699, 41455, 57919,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 70807, 3, 42967,
                                                                       42995, 58927, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 70852, 3, 42995,
                                                                       43023, 58963, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 70897, 3, 43023,
                                                                       43051, 58999, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 70942, 3, 43051,
                                                                       43079, 59035, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 70987, 3, 43079,
                                                                       43107, 59071, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 71032, 3, 43107,
                                                                       43135, 59107, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 71077, 3, 43135,
                                                                       43163, 59143, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 71122, 3, 43163,
                                                                       43191, 59179, ncols,
                                                                       gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 71167, 0, 3,
                                                                       70807, 58927, 70852,
                                                                       43247, 43331, 59215,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 71302, 0, 3,
                                                                       70852, 58963, 70897,
                                                                       43331, 43415, 59323,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 71437, 0, 3,
                                                                       70897, 58999, 70942,
                                                                       43415, 43499, 59431,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 71572, 0, 3,
                                                                       70942, 59035, 70987,
                                                                       43499, 43583, 59539,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 71707, 0, 3,
                                                                       70987, 59071, 71032,
                                                                       43583, 43667, 59647,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 71842, 0, 3,
                                                                       71032, 59107, 71077,
                                                                       43667, 43751, 59755,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 71977, 0, 3,
                                                                       71077, 59143, 71122,
                                                                       43751, 43835, 59863,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 72112, 0, 3,
                                                                       71167, 59215, 71302,
                                                                       44003, 44171, 59971,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 72382, 0, 3,
                                                                       71302, 59323, 71437,
                                                                       44171, 44339, 60187,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 72652, 0, 3,
                                                                       71437, 59431, 71572,
                                                                       44339, 44507, 60403,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 72922, 0, 3,
                                                                       71572, 59539, 71707,
                                                                       44507, 44675, 60619,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 73192, 0, 3,
                                                                       71707, 59647, 71842,
                                                                       44675, 44843, 60835,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 73462, 0, 3,
                                                                       71842, 59755, 71977,
                                                                       44843, 45011, 61051,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 73732, 0, 3,
                                                                       72112, 59971, 72382,
                                                                       45347, 45627, 61267,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 74182, 0, 3,
                                                                       72382, 60187, 72652,
                                                                       45627, 45907, 61627,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 74632, 0, 3,
                                                                       72652, 60403, 72922,
                                                                       45907, 46187, 61987,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 75082, 0, 3,
                                                                       72922, 60619, 73192,
                                                                       46187, 46467, 62347,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 75532, 0, 3,
                                                                       73192, 60835, 73462,
                                                                       46467, 46747, 62707,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 75982, 0, 3,
                                                                       73732, 61267, 74182,
                                                                       47307, 47727, 63067,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 76657, 0, 3,
                                                                       74182, 61627, 74632,
                                                                       47727, 48147, 63607,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 77332, 0, 3,
                                                                       74632, 61987, 75082,
                                                                       48147, 48567, 64147,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 78007, 0, 3,
                                                                       75082, 62347, 75532,
                                                                       48567, 48987, 64687,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 78682, 0, 3,
                                                                       75982, 63067, 76657,
                                                                       49827, 50415, 65227,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 79627, 0, 3,
                                                                       76657, 63607, 77332,
                                                                       50415, 51003, 65983,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 80572, 0, 3,
                                                                       77332, 64147, 78007,
                                                                       51003, 51591, 66739,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 81517, 0, 3,
                                                                       78682, 65227, 79627,
                                                                       52767, 53551, 67495,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 82777, 0, 3,
                                                                       79627, 65983, 80572,
                                                                       53551, 54335, 68503,
                                                                       ncols, gamma, p, q);

                    compute_prim_skl_three_center_electron_repulsion_0(buffer, 84037, 0, 3,
                                                                       81517, 67495, 82777,
                                                                       55903, 56911, 69511,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 85657, 81517, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 87393, 84037, 1620, ncols);
                }
            }
        }

        simdtrf::transform_l_inner(buffer, 86917, 85657, 28, 1, nmax);

        simdtrf::transform_l_inner(buffer, 89013, 87393, 36, 1, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 89625, 86917, 89013, 17, nmax);

        simdtrf::transform_i_inner(buffer, 91053, 89625, 3, 17, nmax);

        simdtrf::transform_p_outer(values + n * npairs, nvalues, buffer, 91053, 221, nmax);
    }

    for (size_t m = 0; m < 663; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
