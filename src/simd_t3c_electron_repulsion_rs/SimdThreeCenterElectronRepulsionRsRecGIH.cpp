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


#include "SimdThreeCenterElectronRepulsionRsRecGIH.hpp"

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
compute_rs_gih_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_gih_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 204937, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 2574 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 204937, 139540, 13994, dimensions);

        for (size_t i = 0; i < nprim_a; i++)
        {
            for (size_t j = 0; j < nprim_b; j++)
            {
                const auto p = a_exps[i] + b_exps[j];

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
                                                            4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
                                                            15}, ncols, fj, i * nprim_b + j, fq,
                                                            omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 22, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15},
                                                        ncols, fj, i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 38, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 41, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 44, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 47, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 50, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 53, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 56, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 59, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 62, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 65, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 68, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 71, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 74, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 77, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 80, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 83, 0, 3, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 86, 0, 3, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 89, 0, 3, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 92, 0, 3, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 95, 0, 3, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 98, 0, 3, 29, 30,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 101, 0, 3, 30, 31,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 104, 0, 3, 31, 32,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 107, 0, 3, 32, 33,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 110, 0, 3, 33, 34,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 113, 0, 3, 34, 35,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 116, 0, 3, 35, 36,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 119, 0, 3, 36, 37,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 122, 0, 3, 7, 8,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 128, 0, 3, 8, 9,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 134, 0, 3, 9, 10,
                                                                       44, 47, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 140, 0, 3, 10, 11,
                                                                       47, 50, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 146, 0, 3, 11, 12,
                                                                       50, 53, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 152, 0, 3, 12, 13,
                                                                       53, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 158, 0, 3, 13, 14,
                                                                       56, 59, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 164, 0, 3, 14, 15,
                                                                       59, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 170, 0, 3, 15, 16,
                                                                       62, 65, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 176, 0, 3, 16, 17,
                                                                       65, 68, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 182, 0, 3, 17, 18,
                                                                       68, 71, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 188, 0, 3, 18, 19,
                                                                       71, 74, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 194, 0, 3, 19, 20,
                                                                       74, 77, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 200, 0, 3, 23, 24,
                                                                       80, 83, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 206, 0, 3, 24, 25,
                                                                       83, 86, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 212, 0, 3, 25, 26,
                                                                       86, 89, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 218, 0, 3, 26, 27,
                                                                       89, 92, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 224, 0, 3, 27, 28,
                                                                       92, 95, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 230, 0, 3, 28, 29,
                                                                       95, 98, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 236, 0, 3, 29, 30,
                                                                       98, 101, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 242, 0, 3, 30, 31,
                                                                       101, 104, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 248, 0, 3, 31, 32,
                                                                       104, 107, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 254, 0, 3, 32, 33,
                                                                       107, 110, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 260, 0, 3, 33, 34,
                                                                       110, 113, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 266, 0, 3, 34, 35,
                                                                       113, 116, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 272, 0, 3, 35, 36,
                                                                       116, 119, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 278, 0, 3, 38, 41,
                                                                       122, 128, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 288, 0, 3, 41, 44,
                                                                       128, 134, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 298, 0, 3, 44, 47,
                                                                       134, 140, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 308, 0, 3, 47, 50,
                                                                       140, 146, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 318, 0, 3, 50, 53,
                                                                       146, 152, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 328, 0, 3, 53, 56,
                                                                       152, 158, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 338, 0, 3, 56, 59,
                                                                       158, 164, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 348, 0, 3, 59, 62,
                                                                       164, 170, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 358, 0, 3, 62, 65,
                                                                       170, 176, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 368, 0, 3, 65, 68,
                                                                       176, 182, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 378, 0, 3, 68, 71,
                                                                       182, 188, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 388, 0, 3, 71, 74,
                                                                       188, 194, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 398, 0, 3, 80, 83,
                                                                       200, 206, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 408, 0, 3, 83, 86,
                                                                       206, 212, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 418, 0, 3, 86, 89,
                                                                       212, 218, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 428, 0, 3, 89, 92,
                                                                       218, 224, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 438, 0, 3, 92, 95,
                                                                       224, 230, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 448, 0, 3, 95, 98,
                                                                       230, 236, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 458, 0, 3, 98,
                                                                       101, 236, 242, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 468, 0, 3, 101,
                                                                       104, 242, 248, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 478, 0, 3, 104,
                                                                       107, 248, 254, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 488, 0, 3, 107,
                                                                       110, 254, 260, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 498, 0, 3, 110,
                                                                       113, 260, 266, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 508, 0, 3, 113,
                                                                       116, 266, 272, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 518, 0, 3, 122,
                                                                       128, 278, 288, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 533, 0, 3, 128,
                                                                       134, 288, 298, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 548, 0, 3, 134,
                                                                       140, 298, 308, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 563, 0, 3, 140,
                                                                       146, 308, 318, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 578, 0, 3, 146,
                                                                       152, 318, 328, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 593, 0, 3, 152,
                                                                       158, 328, 338, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 608, 0, 3, 158,
                                                                       164, 338, 348, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 623, 0, 3, 164,
                                                                       170, 348, 358, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 638, 0, 3, 170,
                                                                       176, 358, 368, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 653, 0, 3, 176,
                                                                       182, 368, 378, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 668, 0, 3, 182,
                                                                       188, 378, 388, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 683, 0, 3, 200,
                                                                       206, 398, 408, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 698, 0, 3, 206,
                                                                       212, 408, 418, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 713, 0, 3, 212,
                                                                       218, 418, 428, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 728, 0, 3, 218,
                                                                       224, 428, 438, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 743, 0, 3, 224,
                                                                       230, 438, 448, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 758, 0, 3, 230,
                                                                       236, 448, 458, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 773, 0, 3, 236,
                                                                       242, 458, 468, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 788, 0, 3, 242,
                                                                       248, 468, 478, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 803, 0, 3, 248,
                                                                       254, 478, 488, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 818, 0, 3, 254,
                                                                       260, 488, 498, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 833, 0, 3, 260,
                                                                       266, 498, 508, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 848, 0, 3, 278,
                                                                       288, 518, 533, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 869, 0, 3, 288,
                                                                       298, 533, 548, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 890, 0, 3, 298,
                                                                       308, 548, 563, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 911, 0, 3, 308,
                                                                       318, 563, 578, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 932, 0, 3, 318,
                                                                       328, 578, 593, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 953, 0, 3, 328,
                                                                       338, 593, 608, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 974, 0, 3, 338,
                                                                       348, 608, 623, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 995, 0, 3, 348,
                                                                       358, 623, 638, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1016, 0, 3, 358,
                                                                       368, 638, 653, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1037, 0, 3, 368,
                                                                       378, 653, 668, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1058, 0, 3, 398,
                                                                       408, 683, 698, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1079, 0, 3, 408,
                                                                       418, 698, 713, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1100, 0, 3, 418,
                                                                       428, 713, 728, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1121, 0, 3, 428,
                                                                       438, 728, 743, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1142, 0, 3, 438,
                                                                       448, 743, 758, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1163, 0, 3, 448,
                                                                       458, 758, 773, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1184, 0, 3, 458,
                                                                       468, 773, 788, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1205, 0, 3, 468,
                                                                       478, 788, 803, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1226, 0, 3, 478,
                                                                       488, 803, 818, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1247, 0, 3, 488,
                                                                       498, 818, 833, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1268, 0, 3, 518,
                                                                       533, 848, 869, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1296, 0, 3, 533,
                                                                       548, 869, 890, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1324, 0, 3, 548,
                                                                       563, 890, 911, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1352, 0, 3, 563,
                                                                       578, 911, 932, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1380, 0, 3, 578,
                                                                       593, 932, 953, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1408, 0, 3, 593,
                                                                       608, 953, 974, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1436, 0, 3, 608,
                                                                       623, 974, 995, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1464, 0, 3, 623,
                                                                       638, 995, 1016, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1492, 0, 3, 638,
                                                                       653, 1016, 1037, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1520, 0, 3, 683,
                                                                       698, 1058, 1079, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1548, 0, 3, 698,
                                                                       713, 1079, 1100, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1576, 0, 3, 713,
                                                                       728, 1100, 1121, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1604, 0, 3, 728,
                                                                       743, 1121, 1142, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1632, 0, 3, 743,
                                                                       758, 1142, 1163, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1660, 0, 3, 758,
                                                                       773, 1163, 1184, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1688, 0, 3, 773,
                                                                       788, 1184, 1205, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1716, 0, 3, 788,
                                                                       803, 1205, 1226, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1744, 0, 3, 803,
                                                                       818, 1226, 1247, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1772, 0, 3, 848,
                                                                       869, 1268, 1296, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1808, 0, 3, 869,
                                                                       890, 1296, 1324, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1844, 0, 3, 890,
                                                                       911, 1324, 1352, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1880, 0, 3, 911,
                                                                       932, 1352, 1380, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1916, 0, 3, 932,
                                                                       953, 1380, 1408, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1952, 0, 3, 953,
                                                                       974, 1408, 1436, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1988, 0, 3, 974,
                                                                       995, 1436, 1464, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2024, 0, 3, 995,
                                                                       1016, 1464, 1492, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2060, 0, 3, 1058,
                                                                       1079, 1520, 1548, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2096, 0, 3, 1079,
                                                                       1100, 1548, 1576, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2132, 0, 3, 1100,
                                                                       1121, 1576, 1604, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2168, 0, 3, 1121,
                                                                       1142, 1604, 1632, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2204, 0, 3, 1142,
                                                                       1163, 1632, 1660, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2240, 0, 3, 1163,
                                                                       1184, 1660, 1688, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2276, 0, 3, 1184,
                                                                       1205, 1688, 1716, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2312, 0, 3, 1205,
                                                                       1226, 1716, 1744, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2348, 0, 3, 1268,
                                                                       1296, 1772, 1808, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2393, 0, 3, 1296,
                                                                       1324, 1808, 1844, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2438, 0, 3, 1324,
                                                                       1352, 1844, 1880, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2483, 0, 3, 1352,
                                                                       1380, 1880, 1916, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2528, 0, 3, 1380,
                                                                       1408, 1916, 1952, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2573, 0, 3, 1408,
                                                                       1436, 1952, 1988, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2618, 0, 3, 1436,
                                                                       1464, 1988, 2024, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2663, 0, 3, 1520,
                                                                       1548, 2060, 2096, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2708, 0, 3, 1548,
                                                                       1576, 2096, 2132, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2753, 0, 3, 1576,
                                                                       1604, 2132, 2168, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2798, 0, 3, 1604,
                                                                       1632, 2168, 2204, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2843, 0, 3, 1632,
                                                                       1660, 2204, 2240, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2888, 0, 3, 1660,
                                                                       1688, 2240, 2276, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2933, 0, 3, 1688,
                                                                       1716, 2276, 2312, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2978, 0, 3, 1772,
                                                                       1808, 2348, 2393, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3033, 0, 3, 1808,
                                                                       1844, 2393, 2438, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3088, 0, 3, 1844,
                                                                       1880, 2438, 2483, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3143, 0, 3, 1880,
                                                                       1916, 2483, 2528, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3198, 0, 3, 1916,
                                                                       1952, 2528, 2573, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3253, 0, 3, 1952,
                                                                       1988, 2573, 2618, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3308, 0, 3, 2060,
                                                                       2096, 2663, 2708, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3363, 0, 3, 2096,
                                                                       2132, 2708, 2753, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3418, 0, 3, 2132,
                                                                       2168, 2753, 2798, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3473, 0, 3, 2168,
                                                                       2204, 2798, 2843, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3528, 0, 3, 2204,
                                                                       2240, 2843, 2888, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3583, 0, 3, 2240,
                                                                       2276, 2888, 2933, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 3638, 0, 3, 2348,
                                                                       2393, 2978, 3033, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 3704, 0, 3, 2393,
                                                                       2438, 3033, 3088, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 3770, 0, 3, 2438,
                                                                       2483, 3088, 3143, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 3836, 0, 3, 2483,
                                                                       2528, 3143, 3198, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 3902, 0, 3, 2528,
                                                                       2573, 3198, 3253, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 3968, 0, 3, 2663,
                                                                       2708, 3308, 3363, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 4034, 0, 3, 2708,
                                                                       2753, 3363, 3418, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 4100, 0, 3, 2753,
                                                                       2798, 3418, 3473, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 4166, 0, 3, 2798,
                                                                       2843, 3473, 3528, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 4232, 0, 3, 2843,
                                                                       2888, 3528, 3583, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4298, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4301, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4304, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4307, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4310, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4313, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4316, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4319, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4322, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4325, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4328, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4331, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4334, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4337, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4340, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4343, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4346, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4349, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4352, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4355, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4358, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4361, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4364, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4367, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4370, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4373, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4376, 3, 34,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4379, 3, 35,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4382, 3, 36,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4385, 3, 37,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4388, 3, 9, 44,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4397, 3, 10, 47,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4406, 3, 11, 50,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4415, 3, 12, 53,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4424, 3, 13, 56,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4433, 3, 14, 59,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4442, 3, 15, 62,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4451, 3, 16, 65,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4460, 3, 17, 68,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4469, 3, 18, 71,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4478, 3, 19, 74,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4487, 3, 20, 77,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4496, 3, 25, 86,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4505, 3, 26, 89,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4514, 3, 27, 92,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4523, 3, 28, 95,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4532, 3, 29, 98,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4541, 3, 30, 101,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4550, 3, 31, 104,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4559, 3, 32, 107,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4568, 3, 33, 110,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4577, 3, 34, 113,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4586, 3, 35, 116,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4595, 3, 36, 119,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4604, 3, 38, 122,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4622, 3, 41, 128,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4640, 3, 44, 134,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4658, 3, 47, 140,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4676, 3, 50, 146,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4694, 3, 53, 152,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4712, 3, 56, 158,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4730, 3, 59, 164,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4748, 3, 62, 170,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4766, 3, 65, 176,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4784, 3, 68, 182,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4802, 3, 71, 188,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4820, 3, 74, 194,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4838, 3, 80, 200,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4856, 3, 83, 206,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4874, 3, 86, 212,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4892, 3, 89, 218,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4910, 3, 92, 224,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4928, 3, 95, 230,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4946, 3, 98, 236,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4964, 3, 101, 242,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4982, 3, 104, 248,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 5000, 3, 107, 254,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 5018, 3, 110, 260,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 5036, 3, 113, 266,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 5054, 3, 116, 272,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5072, 3, 122, 278,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5102, 3, 128, 288,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5132, 3, 134, 298,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5162, 3, 140, 308,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5192, 3, 146, 318,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5222, 3, 152, 328,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5252, 3, 158, 338,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5282, 3, 164, 348,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5312, 3, 170, 358,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5342, 3, 176, 368,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5372, 3, 182, 378,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5402, 3, 188, 388,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5432, 3, 200, 398,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5462, 3, 206, 408,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5492, 3, 212, 418,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5522, 3, 218, 428,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5552, 3, 224, 438,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5582, 3, 230, 448,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5612, 3, 236, 458,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5642, 3, 242, 468,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5672, 3, 248, 478,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5702, 3, 254, 488,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5732, 3, 260, 498,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5762, 3, 266, 508,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5792, 3, 278, 518,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5837, 3, 288, 533,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5882, 3, 298, 548,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5927, 3, 308, 563,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5972, 3, 318, 578,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6017, 3, 328, 593,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6062, 3, 338, 608,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6107, 3, 348, 623,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6152, 3, 358, 638,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6197, 3, 368, 653,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6242, 3, 378, 668,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6287, 3, 398, 683,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6332, 3, 408, 698,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6377, 3, 418, 713,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6422, 3, 428, 728,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6467, 3, 438, 743,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6512, 3, 448, 758,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6557, 3, 458, 773,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6602, 3, 468, 788,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6647, 3, 478, 803,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6692, 3, 488, 818,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6737, 3, 498, 833,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6782, 3, 518, 848,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6845, 3, 533, 869,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6908, 3, 548, 890,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6971, 3, 563, 911,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7034, 3, 578, 932,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7097, 3, 593, 953,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7160, 3, 608, 974,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7223, 3, 623, 995,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7286, 3, 638,
                                                                       1016, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7349, 3, 653,
                                                                       1037, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7412, 3, 683,
                                                                       1058, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7475, 3, 698,
                                                                       1079, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7538, 3, 713,
                                                                       1100, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7601, 3, 728,
                                                                       1121, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7664, 3, 743,
                                                                       1142, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7727, 3, 758,
                                                                       1163, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7790, 3, 773,
                                                                       1184, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7853, 3, 788,
                                                                       1205, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7916, 3, 803,
                                                                       1226, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7979, 3, 818,
                                                                       1247, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8042, 3, 848,
                                                                       1268, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8126, 3, 869,
                                                                       1296, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8210, 3, 890,
                                                                       1324, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8294, 3, 911,
                                                                       1352, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8378, 3, 932,
                                                                       1380, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8462, 3, 953,
                                                                       1408, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8546, 3, 974,
                                                                       1436, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8630, 3, 995,
                                                                       1464, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8714, 3, 1016,
                                                                       1492, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8798, 3, 1058,
                                                                       1520, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8882, 3, 1079,
                                                                       1548, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8966, 3, 1100,
                                                                       1576, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 9050, 3, 1121,
                                                                       1604, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 9134, 3, 1142,
                                                                       1632, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 9218, 3, 1163,
                                                                       1660, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 9302, 3, 1184,
                                                                       1688, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 9386, 3, 1205,
                                                                       1716, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 9470, 3, 1226,
                                                                       1744, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 9554, 3, 1268,
                                                                       1772, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 9662, 3, 1296,
                                                                       1808, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 9770, 3, 1324,
                                                                       1844, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 9878, 3, 1352,
                                                                       1880, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 9986, 3, 1380,
                                                                       1916, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 10094, 3, 1408,
                                                                       1952, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 10202, 3, 1436,
                                                                       1988, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 10310, 3, 1464,
                                                                       2024, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 10418, 3, 1520,
                                                                       2060, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 10526, 3, 1548,
                                                                       2096, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 10634, 3, 1576,
                                                                       2132, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 10742, 3, 1604,
                                                                       2168, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 10850, 3, 1632,
                                                                       2204, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 10958, 3, 1660,
                                                                       2240, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 11066, 3, 1688,
                                                                       2276, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 11174, 3, 1716,
                                                                       2312, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 11282, 3, 1772,
                                                                       2348, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 11417, 3, 1808,
                                                                       2393, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 11552, 3, 1844,
                                                                       2438, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 11687, 3, 1880,
                                                                       2483, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 11822, 3, 1916,
                                                                       2528, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 11957, 3, 1952,
                                                                       2573, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 12092, 3, 1988,
                                                                       2618, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 12227, 3, 2060,
                                                                       2663, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 12362, 3, 2096,
                                                                       2708, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 12497, 3, 2132,
                                                                       2753, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 12632, 3, 2168,
                                                                       2798, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 12767, 3, 2204,
                                                                       2843, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 12902, 3, 2240,
                                                                       2888, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 13037, 3, 2276,
                                                                       2933, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 13172, 3, 2348,
                                                                       2978, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 13337, 3, 2393,
                                                                       3033, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 13502, 3, 2438,
                                                                       3088, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 13667, 3, 2483,
                                                                       3143, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 13832, 3, 2528,
                                                                       3198, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 13997, 3, 2573,
                                                                       3253, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 14162, 3, 2663,
                                                                       3308, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 14327, 3, 2708,
                                                                       3363, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 14492, 3, 2753,
                                                                       3418, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 14657, 3, 2798,
                                                                       3473, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 14822, 3, 2843,
                                                                       3528, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 14987, 3, 2888,
                                                                       3583, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 15152, 3, 2978,
                                                                       3638, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 15350, 3, 3033,
                                                                       3704, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 15548, 3, 3088,
                                                                       3770, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 15746, 3, 3143,
                                                                       3836, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 15944, 3, 3198,
                                                                       3902, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 16142, 3, 3308,
                                                                       3968, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 16340, 3, 3363,
                                                                       4034, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 16538, 3, 3418,
                                                                       4100, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 16736, 3, 3473,
                                                                       4166, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 16934, 3, 3528,
                                                                       4232, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17132, 3, 7, 8,
                                                                       4304, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17138, 3, 8, 9,
                                                                       4307, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17144, 3, 9, 10,
                                                                       4310, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17150, 3, 10, 11,
                                                                       4313, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17156, 3, 11, 12,
                                                                       4316, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17162, 3, 12, 13,
                                                                       4319, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17168, 3, 13, 14,
                                                                       4322, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17174, 3, 14, 15,
                                                                       4325, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17180, 3, 15, 16,
                                                                       4328, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17186, 3, 16, 17,
                                                                       4331, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17192, 3, 17, 18,
                                                                       4334, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17198, 3, 18, 19,
                                                                       4337, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17204, 3, 19, 20,
                                                                       4340, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17210, 3, 23, 24,
                                                                       4349, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17216, 3, 24, 25,
                                                                       4352, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17222, 3, 25, 26,
                                                                       4355, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17228, 3, 26, 27,
                                                                       4358, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17234, 3, 27, 28,
                                                                       4361, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17240, 3, 28, 29,
                                                                       4364, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17246, 3, 29, 30,
                                                                       4367, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17252, 3, 30, 31,
                                                                       4370, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17258, 3, 31, 32,
                                                                       4373, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17264, 3, 32, 33,
                                                                       4376, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17270, 3, 33, 34,
                                                                       4379, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17276, 3, 34, 35,
                                                                       4382, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17282, 3, 35, 36,
                                                                       4385, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 17288, 0, 3,
                                                                       17132, 4304, 17138, 4388,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 17306, 0, 3,
                                                                       17138, 4307, 17144, 4397,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 17324, 0, 3,
                                                                       17144, 4310, 17150, 4406,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 17342, 0, 3,
                                                                       17150, 4313, 17156, 4415,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 17360, 0, 3,
                                                                       17156, 4316, 17162, 4424,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 17378, 0, 3,
                                                                       17162, 4319, 17168, 4433,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 17396, 0, 3,
                                                                       17168, 4322, 17174, 4442,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 17414, 0, 3,
                                                                       17174, 4325, 17180, 4451,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 17432, 0, 3,
                                                                       17180, 4328, 17186, 4460,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 17450, 0, 3,
                                                                       17186, 4331, 17192, 4469,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 17468, 0, 3,
                                                                       17192, 4334, 17198, 4478,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 17486, 0, 3,
                                                                       17198, 4337, 17204, 4487,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 17504, 0, 3,
                                                                       17210, 4349, 17216, 4496,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 17522, 0, 3,
                                                                       17216, 4352, 17222, 4505,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 17540, 0, 3,
                                                                       17222, 4355, 17228, 4514,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 17558, 0, 3,
                                                                       17228, 4358, 17234, 4523,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 17576, 0, 3,
                                                                       17234, 4361, 17240, 4532,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 17594, 0, 3,
                                                                       17240, 4364, 17246, 4541,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 17612, 0, 3,
                                                                       17246, 4367, 17252, 4550,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 17630, 0, 3,
                                                                       17252, 4370, 17258, 4559,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 17648, 0, 3,
                                                                       17258, 4373, 17264, 4568,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 17666, 0, 3,
                                                                       17264, 4376, 17270, 4577,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 17684, 0, 3,
                                                                       17270, 4379, 17276, 4586,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 17702, 0, 3,
                                                                       17276, 4382, 17282, 4595,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17720, 0, 3,
                                                                       17288, 4388, 17306, 122,
                                                                       128, 4640, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17756, 0, 3,
                                                                       17306, 4397, 17324, 128,
                                                                       134, 4658, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17792, 0, 3,
                                                                       17324, 4406, 17342, 134,
                                                                       140, 4676, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17828, 0, 3,
                                                                       17342, 4415, 17360, 140,
                                                                       146, 4694, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17864, 0, 3,
                                                                       17360, 4424, 17378, 146,
                                                                       152, 4712, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17900, 0, 3,
                                                                       17378, 4433, 17396, 152,
                                                                       158, 4730, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17936, 0, 3,
                                                                       17396, 4442, 17414, 158,
                                                                       164, 4748, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17972, 0, 3,
                                                                       17414, 4451, 17432, 164,
                                                                       170, 4766, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 18008, 0, 3,
                                                                       17432, 4460, 17450, 170,
                                                                       176, 4784, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 18044, 0, 3,
                                                                       17450, 4469, 17468, 176,
                                                                       182, 4802, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 18080, 0, 3,
                                                                       17468, 4478, 17486, 182,
                                                                       188, 4820, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 18116, 0, 3,
                                                                       17504, 4496, 17522, 200,
                                                                       206, 4874, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 18152, 0, 3,
                                                                       17522, 4505, 17540, 206,
                                                                       212, 4892, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 18188, 0, 3,
                                                                       17540, 4514, 17558, 212,
                                                                       218, 4910, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 18224, 0, 3,
                                                                       17558, 4523, 17576, 218,
                                                                       224, 4928, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 18260, 0, 3,
                                                                       17576, 4532, 17594, 224,
                                                                       230, 4946, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 18296, 0, 3,
                                                                       17594, 4541, 17612, 230,
                                                                       236, 4964, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 18332, 0, 3,
                                                                       17612, 4550, 17630, 236,
                                                                       242, 4982, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 18368, 0, 3,
                                                                       17630, 4559, 17648, 242,
                                                                       248, 5000, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 18404, 0, 3,
                                                                       17648, 4568, 17666, 248,
                                                                       254, 5018, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 18440, 0, 3,
                                                                       17666, 4577, 17684, 254,
                                                                       260, 5036, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 18476, 0, 3,
                                                                       17684, 4586, 17702, 260,
                                                                       266, 5054, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 18512, 0, 3,
                                                                       17720, 4640, 17756, 278,
                                                                       288, 5132, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 18572, 0, 3,
                                                                       17756, 4658, 17792, 288,
                                                                       298, 5162, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 18632, 0, 3,
                                                                       17792, 4676, 17828, 298,
                                                                       308, 5192, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 18692, 0, 3,
                                                                       17828, 4694, 17864, 308,
                                                                       318, 5222, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 18752, 0, 3,
                                                                       17864, 4712, 17900, 318,
                                                                       328, 5252, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 18812, 0, 3,
                                                                       17900, 4730, 17936, 328,
                                                                       338, 5282, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 18872, 0, 3,
                                                                       17936, 4748, 17972, 338,
                                                                       348, 5312, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 18932, 0, 3,
                                                                       17972, 4766, 18008, 348,
                                                                       358, 5342, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 18992, 0, 3,
                                                                       18008, 4784, 18044, 358,
                                                                       368, 5372, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 19052, 0, 3,
                                                                       18044, 4802, 18080, 368,
                                                                       378, 5402, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 19112, 0, 3,
                                                                       18116, 4874, 18152, 398,
                                                                       408, 5492, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 19172, 0, 3,
                                                                       18152, 4892, 18188, 408,
                                                                       418, 5522, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 19232, 0, 3,
                                                                       18188, 4910, 18224, 418,
                                                                       428, 5552, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 19292, 0, 3,
                                                                       18224, 4928, 18260, 428,
                                                                       438, 5582, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 19352, 0, 3,
                                                                       18260, 4946, 18296, 438,
                                                                       448, 5612, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 19412, 0, 3,
                                                                       18296, 4964, 18332, 448,
                                                                       458, 5642, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 19472, 0, 3,
                                                                       18332, 4982, 18368, 458,
                                                                       468, 5672, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 19532, 0, 3,
                                                                       18368, 5000, 18404, 468,
                                                                       478, 5702, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 19592, 0, 3,
                                                                       18404, 5018, 18440, 478,
                                                                       488, 5732, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 19652, 0, 3,
                                                                       18440, 5036, 18476, 488,
                                                                       498, 5762, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 19712, 0, 3,
                                                                       18512, 5132, 18572, 518,
                                                                       533, 5882, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 19802, 0, 3,
                                                                       18572, 5162, 18632, 533,
                                                                       548, 5927, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 19892, 0, 3,
                                                                       18632, 5192, 18692, 548,
                                                                       563, 5972, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 19982, 0, 3,
                                                                       18692, 5222, 18752, 563,
                                                                       578, 6017, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 20072, 0, 3,
                                                                       18752, 5252, 18812, 578,
                                                                       593, 6062, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 20162, 0, 3,
                                                                       18812, 5282, 18872, 593,
                                                                       608, 6107, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 20252, 0, 3,
                                                                       18872, 5312, 18932, 608,
                                                                       623, 6152, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 20342, 0, 3,
                                                                       18932, 5342, 18992, 623,
                                                                       638, 6197, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 20432, 0, 3,
                                                                       18992, 5372, 19052, 638,
                                                                       653, 6242, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 20522, 0, 3,
                                                                       19112, 5492, 19172, 683,
                                                                       698, 6377, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 20612, 0, 3,
                                                                       19172, 5522, 19232, 698,
                                                                       713, 6422, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 20702, 0, 3,
                                                                       19232, 5552, 19292, 713,
                                                                       728, 6467, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 20792, 0, 3,
                                                                       19292, 5582, 19352, 728,
                                                                       743, 6512, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 20882, 0, 3,
                                                                       19352, 5612, 19412, 743,
                                                                       758, 6557, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 20972, 0, 3,
                                                                       19412, 5642, 19472, 758,
                                                                       773, 6602, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 21062, 0, 3,
                                                                       19472, 5672, 19532, 773,
                                                                       788, 6647, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 21152, 0, 3,
                                                                       19532, 5702, 19592, 788,
                                                                       803, 6692, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 21242, 0, 3,
                                                                       19592, 5732, 19652, 803,
                                                                       818, 6737, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 21332, 0, 3,
                                                                       19712, 5882, 19802, 848,
                                                                       869, 6908, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 21458, 0, 3,
                                                                       19802, 5927, 19892, 869,
                                                                       890, 6971, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 21584, 0, 3,
                                                                       19892, 5972, 19982, 890,
                                                                       911, 7034, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 21710, 0, 3,
                                                                       19982, 6017, 20072, 911,
                                                                       932, 7097, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 21836, 0, 3,
                                                                       20072, 6062, 20162, 932,
                                                                       953, 7160, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 21962, 0, 3,
                                                                       20162, 6107, 20252, 953,
                                                                       974, 7223, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 22088, 0, 3,
                                                                       20252, 6152, 20342, 974,
                                                                       995, 7286, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 22214, 0, 3,
                                                                       20342, 6197, 20432, 995,
                                                                       1016, 7349, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 22340, 0, 3,
                                                                       20522, 6377, 20612, 1058,
                                                                       1079, 7538, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 22466, 0, 3,
                                                                       20612, 6422, 20702, 1079,
                                                                       1100, 7601, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 22592, 0, 3,
                                                                       20702, 6467, 20792, 1100,
                                                                       1121, 7664, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 22718, 0, 3,
                                                                       20792, 6512, 20882, 1121,
                                                                       1142, 7727, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 22844, 0, 3,
                                                                       20882, 6557, 20972, 1142,
                                                                       1163, 7790, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 22970, 0, 3,
                                                                       20972, 6602, 21062, 1163,
                                                                       1184, 7853, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 23096, 0, 3,
                                                                       21062, 6647, 21152, 1184,
                                                                       1205, 7916, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 23222, 0, 3,
                                                                       21152, 6692, 21242, 1205,
                                                                       1226, 7979, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 23348, 0, 3,
                                                                       21332, 6908, 21458, 1268,
                                                                       1296, 8210, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 23516, 0, 3,
                                                                       21458, 6971, 21584, 1296,
                                                                       1324, 8294, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 23684, 0, 3,
                                                                       21584, 7034, 21710, 1324,
                                                                       1352, 8378, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 23852, 0, 3,
                                                                       21710, 7097, 21836, 1352,
                                                                       1380, 8462, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 24020, 0, 3,
                                                                       21836, 7160, 21962, 1380,
                                                                       1408, 8546, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 24188, 0, 3,
                                                                       21962, 7223, 22088, 1408,
                                                                       1436, 8630, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 24356, 0, 3,
                                                                       22088, 7286, 22214, 1436,
                                                                       1464, 8714, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 24524, 0, 3,
                                                                       22340, 7538, 22466, 1520,
                                                                       1548, 8966, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 24692, 0, 3,
                                                                       22466, 7601, 22592, 1548,
                                                                       1576, 9050, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 24860, 0, 3,
                                                                       22592, 7664, 22718, 1576,
                                                                       1604, 9134, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 25028, 0, 3,
                                                                       22718, 7727, 22844, 1604,
                                                                       1632, 9218, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 25196, 0, 3,
                                                                       22844, 7790, 22970, 1632,
                                                                       1660, 9302, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 25364, 0, 3,
                                                                       22970, 7853, 23096, 1660,
                                                                       1688, 9386, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 25532, 0, 3,
                                                                       23096, 7916, 23222, 1688,
                                                                       1716, 9470, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 25700, 0, 3,
                                                                       23348, 8210, 23516, 1772,
                                                                       1808, 9770, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 25916, 0, 3,
                                                                       23516, 8294, 23684, 1808,
                                                                       1844, 9878, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 26132, 0, 3,
                                                                       23684, 8378, 23852, 1844,
                                                                       1880, 9986, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 26348, 0, 3,
                                                                       23852, 8462, 24020, 1880,
                                                                       1916, 10094, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 26564, 0, 3,
                                                                       24020, 8546, 24188, 1916,
                                                                       1952, 10202, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 26780, 0, 3,
                                                                       24188, 8630, 24356, 1952,
                                                                       1988, 10310, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 26996, 0, 3,
                                                                       24524, 8966, 24692, 2060,
                                                                       2096, 10634, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 27212, 0, 3,
                                                                       24692, 9050, 24860, 2096,
                                                                       2132, 10742, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 27428, 0, 3,
                                                                       24860, 9134, 25028, 2132,
                                                                       2168, 10850, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 27644, 0, 3,
                                                                       25028, 9218, 25196, 2168,
                                                                       2204, 10958, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 27860, 0, 3,
                                                                       25196, 9302, 25364, 2204,
                                                                       2240, 11066, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 28076, 0, 3,
                                                                       25364, 9386, 25532, 2240,
                                                                       2276, 11174, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 28292, 0, 3,
                                                                       25700, 9770, 25916, 2348,
                                                                       2393, 11552, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 28562, 0, 3,
                                                                       25916, 9878, 26132, 2393,
                                                                       2438, 11687, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 28832, 0, 3,
                                                                       26132, 9986, 26348, 2438,
                                                                       2483, 11822, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 29102, 0, 3,
                                                                       26348, 10094, 26564, 2483,
                                                                       2528, 11957, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 29372, 0, 3,
                                                                       26564, 10202, 26780, 2528,
                                                                       2573, 12092, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 29642, 0, 3,
                                                                       26996, 10634, 27212, 2663,
                                                                       2708, 12497, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 29912, 0, 3,
                                                                       27212, 10742, 27428, 2708,
                                                                       2753, 12632, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 30182, 0, 3,
                                                                       27428, 10850, 27644, 2753,
                                                                       2798, 12767, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 30452, 0, 3,
                                                                       27644, 10958, 27860, 2798,
                                                                       2843, 12902, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 30722, 0, 3,
                                                                       27860, 11066, 28076, 2843,
                                                                       2888, 13037, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 30992, 0, 3,
                                                                       28292, 11552, 28562, 2978,
                                                                       3033, 13502, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 31322, 0, 3,
                                                                       28562, 11687, 28832, 3033,
                                                                       3088, 13667, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 31652, 0, 3,
                                                                       28832, 11822, 29102, 3088,
                                                                       3143, 13832, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 31982, 0, 3,
                                                                       29102, 11957, 29372, 3143,
                                                                       3198, 13997, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 32312, 0, 3,
                                                                       29642, 12497, 29912, 3308,
                                                                       3363, 14492, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 32642, 0, 3,
                                                                       29912, 12632, 30182, 3363,
                                                                       3418, 14657, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 32972, 0, 3,
                                                                       30182, 12767, 30452, 3418,
                                                                       3473, 14822, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 33302, 0, 3,
                                                                       30452, 12902, 30722, 3473,
                                                                       3528, 14987, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 33632, 0, 3,
                                                                       30992, 13502, 31322, 3638,
                                                                       3704, 15548, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 34028, 0, 3,
                                                                       31322, 13667, 31652, 3704,
                                                                       3770, 15746, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 34424, 0, 3,
                                                                       31652, 13832, 31982, 3770,
                                                                       3836, 15944, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 34820, 0, 3,
                                                                       32312, 14492, 32642, 3968,
                                                                       4034, 16538, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 35216, 0, 3,
                                                                       32642, 14657, 32972, 4034,
                                                                       4100, 16736, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 35612, 0, 3,
                                                                       32972, 14822, 33302, 4100,
                                                                       4166, 16934, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36008, 3, 4298,
                                                                       4301, 17132, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36018, 3, 4301,
                                                                       4304, 17138, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36028, 3, 4304,
                                                                       4307, 17144, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36038, 3, 4307,
                                                                       4310, 17150, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36048, 3, 4310,
                                                                       4313, 17156, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36058, 3, 4313,
                                                                       4316, 17162, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36068, 3, 4316,
                                                                       4319, 17168, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36078, 3, 4319,
                                                                       4322, 17174, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36088, 3, 4322,
                                                                       4325, 17180, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36098, 3, 4325,
                                                                       4328, 17186, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36108, 3, 4328,
                                                                       4331, 17192, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36118, 3, 4331,
                                                                       4334, 17198, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36128, 3, 4334,
                                                                       4337, 17204, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36138, 3, 4343,
                                                                       4346, 17210, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36148, 3, 4346,
                                                                       4349, 17216, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36158, 3, 4349,
                                                                       4352, 17222, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36168, 3, 4352,
                                                                       4355, 17228, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36178, 3, 4355,
                                                                       4358, 17234, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36188, 3, 4358,
                                                                       4361, 17240, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36198, 3, 4361,
                                                                       4364, 17246, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36208, 3, 4364,
                                                                       4367, 17252, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36218, 3, 4367,
                                                                       4370, 17258, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36228, 3, 4370,
                                                                       4373, 17264, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36238, 3, 4373,
                                                                       4376, 17270, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36248, 3, 4376,
                                                                       4379, 17276, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36258, 3, 4379,
                                                                       4382, 17282, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36268, 0, 3,
                                                                       36008, 17132, 36018,
                                                                       17288, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36298, 0, 3,
                                                                       36018, 17138, 36028,
                                                                       17306, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36328, 0, 3,
                                                                       36028, 17144, 36038,
                                                                       17324, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36358, 0, 3,
                                                                       36038, 17150, 36048,
                                                                       17342, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36388, 0, 3,
                                                                       36048, 17156, 36058,
                                                                       17360, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36418, 0, 3,
                                                                       36058, 17162, 36068,
                                                                       17378, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36448, 0, 3,
                                                                       36068, 17168, 36078,
                                                                       17396, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36478, 0, 3,
                                                                       36078, 17174, 36088,
                                                                       17414, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36508, 0, 3,
                                                                       36088, 17180, 36098,
                                                                       17432, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36538, 0, 3,
                                                                       36098, 17186, 36108,
                                                                       17450, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36568, 0, 3,
                                                                       36108, 17192, 36118,
                                                                       17468, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36598, 0, 3,
                                                                       36118, 17198, 36128,
                                                                       17486, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36628, 0, 3,
                                                                       36138, 17210, 36148,
                                                                       17504, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36658, 0, 3,
                                                                       36148, 17216, 36158,
                                                                       17522, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36688, 0, 3,
                                                                       36158, 17222, 36168,
                                                                       17540, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36718, 0, 3,
                                                                       36168, 17228, 36178,
                                                                       17558, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36748, 0, 3,
                                                                       36178, 17234, 36188,
                                                                       17576, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36778, 0, 3,
                                                                       36188, 17240, 36198,
                                                                       17594, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36808, 0, 3,
                                                                       36198, 17246, 36208,
                                                                       17612, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36838, 0, 3,
                                                                       36208, 17252, 36218,
                                                                       17630, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36868, 0, 3,
                                                                       36218, 17258, 36228,
                                                                       17648, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36898, 0, 3,
                                                                       36228, 17264, 36238,
                                                                       17666, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36928, 0, 3,
                                                                       36238, 17270, 36248,
                                                                       17684, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36958, 0, 3,
                                                                       36248, 17276, 36258,
                                                                       17702, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 36988, 0, 3,
                                                                       36268, 17288, 36298, 4604,
                                                                       4622, 17720, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37048, 0, 3,
                                                                       36298, 17306, 36328, 4622,
                                                                       4640, 17756, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37108, 0, 3,
                                                                       36328, 17324, 36358, 4640,
                                                                       4658, 17792, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37168, 0, 3,
                                                                       36358, 17342, 36388, 4658,
                                                                       4676, 17828, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37228, 0, 3,
                                                                       36388, 17360, 36418, 4676,
                                                                       4694, 17864, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37288, 0, 3,
                                                                       36418, 17378, 36448, 4694,
                                                                       4712, 17900, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37348, 0, 3,
                                                                       36448, 17396, 36478, 4712,
                                                                       4730, 17936, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37408, 0, 3,
                                                                       36478, 17414, 36508, 4730,
                                                                       4748, 17972, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37468, 0, 3,
                                                                       36508, 17432, 36538, 4748,
                                                                       4766, 18008, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37528, 0, 3,
                                                                       36538, 17450, 36568, 4766,
                                                                       4784, 18044, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37588, 0, 3,
                                                                       36568, 17468, 36598, 4784,
                                                                       4802, 18080, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37648, 0, 3,
                                                                       36628, 17504, 36658, 4838,
                                                                       4856, 18116, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37708, 0, 3,
                                                                       36658, 17522, 36688, 4856,
                                                                       4874, 18152, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37768, 0, 3,
                                                                       36688, 17540, 36718, 4874,
                                                                       4892, 18188, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37828, 0, 3,
                                                                       36718, 17558, 36748, 4892,
                                                                       4910, 18224, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37888, 0, 3,
                                                                       36748, 17576, 36778, 4910,
                                                                       4928, 18260, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37948, 0, 3,
                                                                       36778, 17594, 36808, 4928,
                                                                       4946, 18296, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 38008, 0, 3,
                                                                       36808, 17612, 36838, 4946,
                                                                       4964, 18332, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 38068, 0, 3,
                                                                       36838, 17630, 36868, 4964,
                                                                       4982, 18368, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 38128, 0, 3,
                                                                       36868, 17648, 36898, 4982,
                                                                       5000, 18404, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 38188, 0, 3,
                                                                       36898, 17666, 36928, 5000,
                                                                       5018, 18440, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 38248, 0, 3,
                                                                       36928, 17684, 36958, 5018,
                                                                       5036, 18476, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 38308, 0, 3,
                                                                       36988, 17720, 37048, 5072,
                                                                       5102, 18512, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 38408, 0, 3,
                                                                       37048, 17756, 37108, 5102,
                                                                       5132, 18572, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 38508, 0, 3,
                                                                       37108, 17792, 37168, 5132,
                                                                       5162, 18632, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 38608, 0, 3,
                                                                       37168, 17828, 37228, 5162,
                                                                       5192, 18692, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 38708, 0, 3,
                                                                       37228, 17864, 37288, 5192,
                                                                       5222, 18752, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 38808, 0, 3,
                                                                       37288, 17900, 37348, 5222,
                                                                       5252, 18812, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 38908, 0, 3,
                                                                       37348, 17936, 37408, 5252,
                                                                       5282, 18872, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 39008, 0, 3,
                                                                       37408, 17972, 37468, 5282,
                                                                       5312, 18932, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 39108, 0, 3,
                                                                       37468, 18008, 37528, 5312,
                                                                       5342, 18992, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 39208, 0, 3,
                                                                       37528, 18044, 37588, 5342,
                                                                       5372, 19052, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 39308, 0, 3,
                                                                       37648, 18116, 37708, 5432,
                                                                       5462, 19112, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 39408, 0, 3,
                                                                       37708, 18152, 37768, 5462,
                                                                       5492, 19172, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 39508, 0, 3,
                                                                       37768, 18188, 37828, 5492,
                                                                       5522, 19232, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 39608, 0, 3,
                                                                       37828, 18224, 37888, 5522,
                                                                       5552, 19292, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 39708, 0, 3,
                                                                       37888, 18260, 37948, 5552,
                                                                       5582, 19352, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 39808, 0, 3,
                                                                       37948, 18296, 38008, 5582,
                                                                       5612, 19412, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 39908, 0, 3,
                                                                       38008, 18332, 38068, 5612,
                                                                       5642, 19472, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 40008, 0, 3,
                                                                       38068, 18368, 38128, 5642,
                                                                       5672, 19532, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 40108, 0, 3,
                                                                       38128, 18404, 38188, 5672,
                                                                       5702, 19592, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 40208, 0, 3,
                                                                       38188, 18440, 38248, 5702,
                                                                       5732, 19652, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 40308, 0, 3,
                                                                       38308, 18512, 38408, 5792,
                                                                       5837, 19712, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 40458, 0, 3,
                                                                       38408, 18572, 38508, 5837,
                                                                       5882, 19802, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 40608, 0, 3,
                                                                       38508, 18632, 38608, 5882,
                                                                       5927, 19892, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 40758, 0, 3,
                                                                       38608, 18692, 38708, 5927,
                                                                       5972, 19982, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 40908, 0, 3,
                                                                       38708, 18752, 38808, 5972,
                                                                       6017, 20072, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 41058, 0, 3,
                                                                       38808, 18812, 38908, 6017,
                                                                       6062, 20162, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 41208, 0, 3,
                                                                       38908, 18872, 39008, 6062,
                                                                       6107, 20252, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 41358, 0, 3,
                                                                       39008, 18932, 39108, 6107,
                                                                       6152, 20342, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 41508, 0, 3,
                                                                       39108, 18992, 39208, 6152,
                                                                       6197, 20432, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 41658, 0, 3,
                                                                       39308, 19112, 39408, 6287,
                                                                       6332, 20522, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 41808, 0, 3,
                                                                       39408, 19172, 39508, 6332,
                                                                       6377, 20612, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 41958, 0, 3,
                                                                       39508, 19232, 39608, 6377,
                                                                       6422, 20702, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 42108, 0, 3,
                                                                       39608, 19292, 39708, 6422,
                                                                       6467, 20792, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 42258, 0, 3,
                                                                       39708, 19352, 39808, 6467,
                                                                       6512, 20882, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 42408, 0, 3,
                                                                       39808, 19412, 39908, 6512,
                                                                       6557, 20972, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 42558, 0, 3,
                                                                       39908, 19472, 40008, 6557,
                                                                       6602, 21062, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 42708, 0, 3,
                                                                       40008, 19532, 40108, 6602,
                                                                       6647, 21152, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 42858, 0, 3,
                                                                       40108, 19592, 40208, 6647,
                                                                       6692, 21242, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 43008, 0, 3,
                                                                       40308, 19712, 40458, 6782,
                                                                       6845, 21332, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 43218, 0, 3,
                                                                       40458, 19802, 40608, 6845,
                                                                       6908, 21458, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 43428, 0, 3,
                                                                       40608, 19892, 40758, 6908,
                                                                       6971, 21584, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 43638, 0, 3,
                                                                       40758, 19982, 40908, 6971,
                                                                       7034, 21710, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 43848, 0, 3,
                                                                       40908, 20072, 41058, 7034,
                                                                       7097, 21836, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 44058, 0, 3,
                                                                       41058, 20162, 41208, 7097,
                                                                       7160, 21962, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 44268, 0, 3,
                                                                       41208, 20252, 41358, 7160,
                                                                       7223, 22088, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 44478, 0, 3,
                                                                       41358, 20342, 41508, 7223,
                                                                       7286, 22214, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 44688, 0, 3,
                                                                       41658, 20522, 41808, 7412,
                                                                       7475, 22340, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 44898, 0, 3,
                                                                       41808, 20612, 41958, 7475,
                                                                       7538, 22466, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 45108, 0, 3,
                                                                       41958, 20702, 42108, 7538,
                                                                       7601, 22592, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 45318, 0, 3,
                                                                       42108, 20792, 42258, 7601,
                                                                       7664, 22718, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 45528, 0, 3,
                                                                       42258, 20882, 42408, 7664,
                                                                       7727, 22844, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 45738, 0, 3,
                                                                       42408, 20972, 42558, 7727,
                                                                       7790, 22970, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 45948, 0, 3,
                                                                       42558, 21062, 42708, 7790,
                                                                       7853, 23096, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 46158, 0, 3,
                                                                       42708, 21152, 42858, 7853,
                                                                       7916, 23222, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 46368, 0, 3,
                                                                       43008, 21332, 43218, 8042,
                                                                       8126, 23348, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 46648, 0, 3,
                                                                       43218, 21458, 43428, 8126,
                                                                       8210, 23516, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 46928, 0, 3,
                                                                       43428, 21584, 43638, 8210,
                                                                       8294, 23684, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 47208, 0, 3,
                                                                       43638, 21710, 43848, 8294,
                                                                       8378, 23852, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 47488, 0, 3,
                                                                       43848, 21836, 44058, 8378,
                                                                       8462, 24020, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 47768, 0, 3,
                                                                       44058, 21962, 44268, 8462,
                                                                       8546, 24188, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 48048, 0, 3,
                                                                       44268, 22088, 44478, 8546,
                                                                       8630, 24356, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 48328, 0, 3,
                                                                       44688, 22340, 44898, 8798,
                                                                       8882, 24524, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 48608, 0, 3,
                                                                       44898, 22466, 45108, 8882,
                                                                       8966, 24692, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 48888, 0, 3,
                                                                       45108, 22592, 45318, 8966,
                                                                       9050, 24860, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 49168, 0, 3,
                                                                       45318, 22718, 45528, 9050,
                                                                       9134, 25028, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 49448, 0, 3,
                                                                       45528, 22844, 45738, 9134,
                                                                       9218, 25196, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 49728, 0, 3,
                                                                       45738, 22970, 45948, 9218,
                                                                       9302, 25364, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 50008, 0, 3,
                                                                       45948, 23096, 46158, 9302,
                                                                       9386, 25532, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 50288, 0, 3,
                                                                       46368, 23348, 46648, 9554,
                                                                       9662, 25700, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 50648, 0, 3,
                                                                       46648, 23516, 46928, 9662,
                                                                       9770, 25916, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 51008, 0, 3,
                                                                       46928, 23684, 47208, 9770,
                                                                       9878, 26132, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 51368, 0, 3,
                                                                       47208, 23852, 47488, 9878,
                                                                       9986, 26348, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 51728, 0, 3,
                                                                       47488, 24020, 47768, 9986,
                                                                       10094, 26564, ncols,
                                                                       gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 52088, 0, 3,
                                                                       47768, 24188, 48048,
                                                                       10094, 10202, 26780,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 52448, 0, 3,
                                                                       48328, 24524, 48608,
                                                                       10418, 10526, 26996,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 52808, 0, 3,
                                                                       48608, 24692, 48888,
                                                                       10526, 10634, 27212,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 53168, 0, 3,
                                                                       48888, 24860, 49168,
                                                                       10634, 10742, 27428,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 53528, 0, 3,
                                                                       49168, 25028, 49448,
                                                                       10742, 10850, 27644,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 53888, 0, 3,
                                                                       49448, 25196, 49728,
                                                                       10850, 10958, 27860,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 54248, 0, 3,
                                                                       49728, 25364, 50008,
                                                                       10958, 11066, 28076,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 54608, 0, 3,
                                                                       50288, 25700, 50648,
                                                                       11282, 11417, 28292,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 55058, 0, 3,
                                                                       50648, 25916, 51008,
                                                                       11417, 11552, 28562,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 55508, 0, 3,
                                                                       51008, 26132, 51368,
                                                                       11552, 11687, 28832,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 55958, 0, 3,
                                                                       51368, 26348, 51728,
                                                                       11687, 11822, 29102,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 56408, 0, 3,
                                                                       51728, 26564, 52088,
                                                                       11822, 11957, 29372,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 56858, 0, 3,
                                                                       52448, 26996, 52808,
                                                                       12227, 12362, 29642,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 57308, 0, 3,
                                                                       52808, 27212, 53168,
                                                                       12362, 12497, 29912,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 57758, 0, 3,
                                                                       53168, 27428, 53528,
                                                                       12497, 12632, 30182,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 58208, 0, 3,
                                                                       53528, 27644, 53888,
                                                                       12632, 12767, 30452,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 58658, 0, 3,
                                                                       53888, 27860, 54248,
                                                                       12767, 12902, 30722,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 59108, 0, 3,
                                                                       54608, 28292, 55058,
                                                                       13172, 13337, 30992,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 59658, 0, 3,
                                                                       55058, 28562, 55508,
                                                                       13337, 13502, 31322,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 60208, 0, 3,
                                                                       55508, 28832, 55958,
                                                                       13502, 13667, 31652,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 60758, 0, 3,
                                                                       55958, 29102, 56408,
                                                                       13667, 13832, 31982,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 61308, 0, 3,
                                                                       56858, 29642, 57308,
                                                                       14162, 14327, 32312,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 61858, 0, 3,
                                                                       57308, 29912, 57758,
                                                                       14327, 14492, 32642,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 62408, 0, 3,
                                                                       57758, 30182, 58208,
                                                                       14492, 14657, 32972,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 62958, 0, 3,
                                                                       58208, 30452, 58658,
                                                                       14657, 14822, 33302,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 63508, 0, 3,
                                                                       59108, 30992, 59658,
                                                                       15152, 15350, 33632,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 64168, 0, 3,
                                                                       59658, 31322, 60208,
                                                                       15350, 15548, 34028,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 64828, 0, 3,
                                                                       60208, 31652, 60758,
                                                                       15548, 15746, 34424,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 65488, 0, 3,
                                                                       61308, 32312, 61858,
                                                                       16142, 16340, 34820,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 66148, 0, 3,
                                                                       61858, 32642, 62408,
                                                                       16340, 16538, 35216,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 66808, 0, 3,
                                                                       62408, 32972, 62958,
                                                                       16538, 16736, 35612,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67468, 3, 17132,
                                                                       17138, 36028, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67483, 3, 17138,
                                                                       17144, 36038, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67498, 3, 17144,
                                                                       17150, 36048, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67513, 3, 17150,
                                                                       17156, 36058, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67528, 3, 17156,
                                                                       17162, 36068, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67543, 3, 17162,
                                                                       17168, 36078, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67558, 3, 17168,
                                                                       17174, 36088, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67573, 3, 17174,
                                                                       17180, 36098, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67588, 3, 17180,
                                                                       17186, 36108, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67603, 3, 17186,
                                                                       17192, 36118, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67618, 3, 17192,
                                                                       17198, 36128, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67633, 3, 17210,
                                                                       17216, 36158, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67648, 3, 17216,
                                                                       17222, 36168, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67663, 3, 17222,
                                                                       17228, 36178, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67678, 3, 17228,
                                                                       17234, 36188, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67693, 3, 17234,
                                                                       17240, 36198, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67708, 3, 17240,
                                                                       17246, 36208, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67723, 3, 17246,
                                                                       17252, 36218, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67738, 3, 17252,
                                                                       17258, 36228, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67753, 3, 17258,
                                                                       17264, 36238, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67768, 3, 17264,
                                                                       17270, 36248, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67783, 3, 17270,
                                                                       17276, 36258, ncols,
                                                                       gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 67798, 0, 3,
                                                                       67468, 36028, 67483,
                                                                       17288, 17306, 36328,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 67843, 0, 3,
                                                                       67483, 36038, 67498,
                                                                       17306, 17324, 36358,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 67888, 0, 3,
                                                                       67498, 36048, 67513,
                                                                       17324, 17342, 36388,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 67933, 0, 3,
                                                                       67513, 36058, 67528,
                                                                       17342, 17360, 36418,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 67978, 0, 3,
                                                                       67528, 36068, 67543,
                                                                       17360, 17378, 36448,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 68023, 0, 3,
                                                                       67543, 36078, 67558,
                                                                       17378, 17396, 36478,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 68068, 0, 3,
                                                                       67558, 36088, 67573,
                                                                       17396, 17414, 36508,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 68113, 0, 3,
                                                                       67573, 36098, 67588,
                                                                       17414, 17432, 36538,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 68158, 0, 3,
                                                                       67588, 36108, 67603,
                                                                       17432, 17450, 36568,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 68203, 0, 3,
                                                                       67603, 36118, 67618,
                                                                       17450, 17468, 36598,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 68248, 0, 3,
                                                                       67633, 36158, 67648,
                                                                       17504, 17522, 36688,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 68293, 0, 3,
                                                                       67648, 36168, 67663,
                                                                       17522, 17540, 36718,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 68338, 0, 3,
                                                                       67663, 36178, 67678,
                                                                       17540, 17558, 36748,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 68383, 0, 3,
                                                                       67678, 36188, 67693,
                                                                       17558, 17576, 36778,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 68428, 0, 3,
                                                                       67693, 36198, 67708,
                                                                       17576, 17594, 36808,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 68473, 0, 3,
                                                                       67708, 36208, 67723,
                                                                       17594, 17612, 36838,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 68518, 0, 3,
                                                                       67723, 36218, 67738,
                                                                       17612, 17630, 36868,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 68563, 0, 3,
                                                                       67738, 36228, 67753,
                                                                       17630, 17648, 36898,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 68608, 0, 3,
                                                                       67753, 36238, 67768,
                                                                       17648, 17666, 36928,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 68653, 0, 3,
                                                                       67768, 36248, 67783,
                                                                       17666, 17684, 36958,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 68698, 0, 3,
                                                                       67798, 36328, 67843,
                                                                       17720, 17756, 37108,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 68788, 0, 3,
                                                                       67843, 36358, 67888,
                                                                       17756, 17792, 37168,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 68878, 0, 3,
                                                                       67888, 36388, 67933,
                                                                       17792, 17828, 37228,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 68968, 0, 3,
                                                                       67933, 36418, 67978,
                                                                       17828, 17864, 37288,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 69058, 0, 3,
                                                                       67978, 36448, 68023,
                                                                       17864, 17900, 37348,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 69148, 0, 3,
                                                                       68023, 36478, 68068,
                                                                       17900, 17936, 37408,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 69238, 0, 3,
                                                                       68068, 36508, 68113,
                                                                       17936, 17972, 37468,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 69328, 0, 3,
                                                                       68113, 36538, 68158,
                                                                       17972, 18008, 37528,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 69418, 0, 3,
                                                                       68158, 36568, 68203,
                                                                       18008, 18044, 37588,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 69508, 0, 3,
                                                                       68248, 36688, 68293,
                                                                       18116, 18152, 37768,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 69598, 0, 3,
                                                                       68293, 36718, 68338,
                                                                       18152, 18188, 37828,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 69688, 0, 3,
                                                                       68338, 36748, 68383,
                                                                       18188, 18224, 37888,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 69778, 0, 3,
                                                                       68383, 36778, 68428,
                                                                       18224, 18260, 37948,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 69868, 0, 3,
                                                                       68428, 36808, 68473,
                                                                       18260, 18296, 38008,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 69958, 0, 3,
                                                                       68473, 36838, 68518,
                                                                       18296, 18332, 38068,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 70048, 0, 3,
                                                                       68518, 36868, 68563,
                                                                       18332, 18368, 38128,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 70138, 0, 3,
                                                                       68563, 36898, 68608,
                                                                       18368, 18404, 38188,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 70228, 0, 3,
                                                                       68608, 36928, 68653,
                                                                       18404, 18440, 38248,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 70318, 0, 3,
                                                                       68698, 37108, 68788,
                                                                       18512, 18572, 38508,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 70468, 0, 3,
                                                                       68788, 37168, 68878,
                                                                       18572, 18632, 38608,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 70618, 0, 3,
                                                                       68878, 37228, 68968,
                                                                       18632, 18692, 38708,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 70768, 0, 3,
                                                                       68968, 37288, 69058,
                                                                       18692, 18752, 38808,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 70918, 0, 3,
                                                                       69058, 37348, 69148,
                                                                       18752, 18812, 38908,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 71068, 0, 3,
                                                                       69148, 37408, 69238,
                                                                       18812, 18872, 39008,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 71218, 0, 3,
                                                                       69238, 37468, 69328,
                                                                       18872, 18932, 39108,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 71368, 0, 3,
                                                                       69328, 37528, 69418,
                                                                       18932, 18992, 39208,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 71518, 0, 3,
                                                                       69508, 37768, 69598,
                                                                       19112, 19172, 39508,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 71668, 0, 3,
                                                                       69598, 37828, 69688,
                                                                       19172, 19232, 39608,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 71818, 0, 3,
                                                                       69688, 37888, 69778,
                                                                       19232, 19292, 39708,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 71968, 0, 3,
                                                                       69778, 37948, 69868,
                                                                       19292, 19352, 39808,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 72118, 0, 3,
                                                                       69868, 38008, 69958,
                                                                       19352, 19412, 39908,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 72268, 0, 3,
                                                                       69958, 38068, 70048,
                                                                       19412, 19472, 40008,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 72418, 0, 3,
                                                                       70048, 38128, 70138,
                                                                       19472, 19532, 40108,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 72568, 0, 3,
                                                                       70138, 38188, 70228,
                                                                       19532, 19592, 40208,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 72718, 0, 3,
                                                                       70318, 38508, 70468,
                                                                       19712, 19802, 40608,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 72943, 0, 3,
                                                                       70468, 38608, 70618,
                                                                       19802, 19892, 40758,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 73168, 0, 3,
                                                                       70618, 38708, 70768,
                                                                       19892, 19982, 40908,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 73393, 0, 3,
                                                                       70768, 38808, 70918,
                                                                       19982, 20072, 41058,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 73618, 0, 3,
                                                                       70918, 38908, 71068,
                                                                       20072, 20162, 41208,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 73843, 0, 3,
                                                                       71068, 39008, 71218,
                                                                       20162, 20252, 41358,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 74068, 0, 3,
                                                                       71218, 39108, 71368,
                                                                       20252, 20342, 41508,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 74293, 0, 3,
                                                                       71518, 39508, 71668,
                                                                       20522, 20612, 41958,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 74518, 0, 3,
                                                                       71668, 39608, 71818,
                                                                       20612, 20702, 42108,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 74743, 0, 3,
                                                                       71818, 39708, 71968,
                                                                       20702, 20792, 42258,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 74968, 0, 3,
                                                                       71968, 39808, 72118,
                                                                       20792, 20882, 42408,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 75193, 0, 3,
                                                                       72118, 39908, 72268,
                                                                       20882, 20972, 42558,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 75418, 0, 3,
                                                                       72268, 40008, 72418,
                                                                       20972, 21062, 42708,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 75643, 0, 3,
                                                                       72418, 40108, 72568,
                                                                       21062, 21152, 42858,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 75868, 0, 3,
                                                                       72718, 40608, 72943,
                                                                       21332, 21458, 43428,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 76183, 0, 3,
                                                                       72943, 40758, 73168,
                                                                       21458, 21584, 43638,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 76498, 0, 3,
                                                                       73168, 40908, 73393,
                                                                       21584, 21710, 43848,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 76813, 0, 3,
                                                                       73393, 41058, 73618,
                                                                       21710, 21836, 44058,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 77128, 0, 3,
                                                                       73618, 41208, 73843,
                                                                       21836, 21962, 44268,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 77443, 0, 3,
                                                                       73843, 41358, 74068,
                                                                       21962, 22088, 44478,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 77758, 0, 3,
                                                                       74293, 41958, 74518,
                                                                       22340, 22466, 45108,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 78073, 0, 3,
                                                                       74518, 42108, 74743,
                                                                       22466, 22592, 45318,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 78388, 0, 3,
                                                                       74743, 42258, 74968,
                                                                       22592, 22718, 45528,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 78703, 0, 3,
                                                                       74968, 42408, 75193,
                                                                       22718, 22844, 45738,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 79018, 0, 3,
                                                                       75193, 42558, 75418,
                                                                       22844, 22970, 45948,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 79333, 0, 3,
                                                                       75418, 42708, 75643,
                                                                       22970, 23096, 46158,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 79648, 0, 3,
                                                                       75868, 43428, 76183,
                                                                       23348, 23516, 46928,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 80068, 0, 3,
                                                                       76183, 43638, 76498,
                                                                       23516, 23684, 47208,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 80488, 0, 3,
                                                                       76498, 43848, 76813,
                                                                       23684, 23852, 47488,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 80908, 0, 3,
                                                                       76813, 44058, 77128,
                                                                       23852, 24020, 47768,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 81328, 0, 3,
                                                                       77128, 44268, 77443,
                                                                       24020, 24188, 48048,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 81748, 0, 3,
                                                                       77758, 45108, 78073,
                                                                       24524, 24692, 48888,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 82168, 0, 3,
                                                                       78073, 45318, 78388,
                                                                       24692, 24860, 49168,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 82588, 0, 3,
                                                                       78388, 45528, 78703,
                                                                       24860, 25028, 49448,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 83008, 0, 3,
                                                                       78703, 45738, 79018,
                                                                       25028, 25196, 49728,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 83428, 0, 3,
                                                                       79018, 45948, 79333,
                                                                       25196, 25364, 50008,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 83848, 0, 3,
                                                                       79648, 46928, 80068,
                                                                       25700, 25916, 51008,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 84388, 0, 3,
                                                                       80068, 47208, 80488,
                                                                       25916, 26132, 51368,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 84928, 0, 3,
                                                                       80488, 47488, 80908,
                                                                       26132, 26348, 51728,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 85468, 0, 3,
                                                                       80908, 47768, 81328,
                                                                       26348, 26564, 52088,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 86008, 0, 3,
                                                                       81748, 48888, 82168,
                                                                       26996, 27212, 53168,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 86548, 0, 3,
                                                                       82168, 49168, 82588,
                                                                       27212, 27428, 53528,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 87088, 0, 3,
                                                                       82588, 49448, 83008,
                                                                       27428, 27644, 53888,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 87628, 0, 3,
                                                                       83008, 49728, 83428,
                                                                       27644, 27860, 54248,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 88168, 0, 3,
                                                                       83848, 51008, 84388,
                                                                       28292, 28562, 55508,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 88843, 0, 3,
                                                                       84388, 51368, 84928,
                                                                       28562, 28832, 55958,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 89518, 0, 3,
                                                                       84928, 51728, 85468,
                                                                       28832, 29102, 56408,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 90193, 0, 3,
                                                                       86008, 53168, 86548,
                                                                       29642, 29912, 57758,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 90868, 0, 3,
                                                                       86548, 53528, 87088,
                                                                       29912, 30182, 58208,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 91543, 0, 3,
                                                                       87088, 53888, 87628,
                                                                       30182, 30452, 58658,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 92218, 0, 3,
                                                                       88168, 55508, 88843,
                                                                       30992, 31322, 60208,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 93043, 0, 3,
                                                                       88843, 55958, 89518,
                                                                       31322, 31652, 60758,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 93868, 0, 3,
                                                                       90193, 57758, 90868,
                                                                       32312, 32642, 62408,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 94693, 0, 3,
                                                                       90868, 58208, 91543,
                                                                       32642, 32972, 62958,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 95518, 0, 3,
                                                                       92218, 60208, 93043,
                                                                       33632, 34028, 64828,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 96508, 0, 3,
                                                                       93868, 62408, 94693,
                                                                       34820, 35216, 66808,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97498, 3, 36008,
                                                                       36018, 67468, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97519, 3, 36018,
                                                                       36028, 67483, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97540, 3, 36028,
                                                                       36038, 67498, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97561, 3, 36038,
                                                                       36048, 67513, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97582, 3, 36048,
                                                                       36058, 67528, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97603, 3, 36058,
                                                                       36068, 67543, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97624, 3, 36068,
                                                                       36078, 67558, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97645, 3, 36078,
                                                                       36088, 67573, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97666, 3, 36088,
                                                                       36098, 67588, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97687, 3, 36098,
                                                                       36108, 67603, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97708, 3, 36108,
                                                                       36118, 67618, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97729, 3, 36138,
                                                                       36148, 67633, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97750, 3, 36148,
                                                                       36158, 67648, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97771, 3, 36158,
                                                                       36168, 67663, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97792, 3, 36168,
                                                                       36178, 67678, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97813, 3, 36178,
                                                                       36188, 67693, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97834, 3, 36188,
                                                                       36198, 67708, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97855, 3, 36198,
                                                                       36208, 67723, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97876, 3, 36208,
                                                                       36218, 67738, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97897, 3, 36218,
                                                                       36228, 67753, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97918, 3, 36228,
                                                                       36238, 67768, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97939, 3, 36238,
                                                                       36248, 67783, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 97960, 0, 3,
                                                                       97498, 67468, 97519,
                                                                       36268, 36298, 67798,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 98023, 0, 3,
                                                                       97519, 67483, 97540,
                                                                       36298, 36328, 67843,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 98086, 0, 3,
                                                                       97540, 67498, 97561,
                                                                       36328, 36358, 67888,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 98149, 0, 3,
                                                                       97561, 67513, 97582,
                                                                       36358, 36388, 67933,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 98212, 0, 3,
                                                                       97582, 67528, 97603,
                                                                       36388, 36418, 67978,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 98275, 0, 3,
                                                                       97603, 67543, 97624,
                                                                       36418, 36448, 68023,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 98338, 0, 3,
                                                                       97624, 67558, 97645,
                                                                       36448, 36478, 68068,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 98401, 0, 3,
                                                                       97645, 67573, 97666,
                                                                       36478, 36508, 68113,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 98464, 0, 3,
                                                                       97666, 67588, 97687,
                                                                       36508, 36538, 68158,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 98527, 0, 3,
                                                                       97687, 67603, 97708,
                                                                       36538, 36568, 68203,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 98590, 0, 3,
                                                                       97729, 67633, 97750,
                                                                       36628, 36658, 68248,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 98653, 0, 3,
                                                                       97750, 67648, 97771,
                                                                       36658, 36688, 68293,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 98716, 0, 3,
                                                                       97771, 67663, 97792,
                                                                       36688, 36718, 68338,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 98779, 0, 3,
                                                                       97792, 67678, 97813,
                                                                       36718, 36748, 68383,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 98842, 0, 3,
                                                                       97813, 67693, 97834,
                                                                       36748, 36778, 68428,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 98905, 0, 3,
                                                                       97834, 67708, 97855,
                                                                       36778, 36808, 68473,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 98968, 0, 3,
                                                                       97855, 67723, 97876,
                                                                       36808, 36838, 68518,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 99031, 0, 3,
                                                                       97876, 67738, 97897,
                                                                       36838, 36868, 68563,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 99094, 0, 3,
                                                                       97897, 67753, 97918,
                                                                       36868, 36898, 68608,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 99157, 0, 3,
                                                                       97918, 67768, 97939,
                                                                       36898, 36928, 68653,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 99220, 0, 3,
                                                                       97960, 67798, 98023,
                                                                       36988, 37048, 68698,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 99346, 0, 3,
                                                                       98023, 67843, 98086,
                                                                       37048, 37108, 68788,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 99472, 0, 3,
                                                                       98086, 67888, 98149,
                                                                       37108, 37168, 68878,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 99598, 0, 3,
                                                                       98149, 67933, 98212,
                                                                       37168, 37228, 68968,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 99724, 0, 3,
                                                                       98212, 67978, 98275,
                                                                       37228, 37288, 69058,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 99850, 0, 3,
                                                                       98275, 68023, 98338,
                                                                       37288, 37348, 69148,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 99976, 0, 3,
                                                                       98338, 68068, 98401,
                                                                       37348, 37408, 69238,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 100102, 0, 3,
                                                                       98401, 68113, 98464,
                                                                       37408, 37468, 69328,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 100228, 0, 3,
                                                                       98464, 68158, 98527,
                                                                       37468, 37528, 69418,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 100354, 0, 3,
                                                                       98590, 68248, 98653,
                                                                       37648, 37708, 69508,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 100480, 0, 3,
                                                                       98653, 68293, 98716,
                                                                       37708, 37768, 69598,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 100606, 0, 3,
                                                                       98716, 68338, 98779,
                                                                       37768, 37828, 69688,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 100732, 0, 3,
                                                                       98779, 68383, 98842,
                                                                       37828, 37888, 69778,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 100858, 0, 3,
                                                                       98842, 68428, 98905,
                                                                       37888, 37948, 69868,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 100984, 0, 3,
                                                                       98905, 68473, 98968,
                                                                       37948, 38008, 69958,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 101110, 0, 3,
                                                                       98968, 68518, 99031,
                                                                       38008, 38068, 70048,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 101236, 0, 3,
                                                                       99031, 68563, 99094,
                                                                       38068, 38128, 70138,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 101362, 0, 3,
                                                                       99094, 68608, 99157,
                                                                       38128, 38188, 70228,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 101488, 0, 3,
                                                                       99220, 68698, 99346,
                                                                       38308, 38408, 70318,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 101698, 0, 3,
                                                                       99346, 68788, 99472,
                                                                       38408, 38508, 70468,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 101908, 0, 3,
                                                                       99472, 68878, 99598,
                                                                       38508, 38608, 70618,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 102118, 0, 3,
                                                                       99598, 68968, 99724,
                                                                       38608, 38708, 70768,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 102328, 0, 3,
                                                                       99724, 69058, 99850,
                                                                       38708, 38808, 70918,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 102538, 0, 3,
                                                                       99850, 69148, 99976,
                                                                       38808, 38908, 71068,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 102748, 0, 3,
                                                                       99976, 69238, 100102,
                                                                       38908, 39008, 71218,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 102958, 0, 3,
                                                                       100102, 69328, 100228,
                                                                       39008, 39108, 71368,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 103168, 0, 3,
                                                                       100354, 69508, 100480,
                                                                       39308, 39408, 71518,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 103378, 0, 3,
                                                                       100480, 69598, 100606,
                                                                       39408, 39508, 71668,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 103588, 0, 3,
                                                                       100606, 69688, 100732,
                                                                       39508, 39608, 71818,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 103798, 0, 3,
                                                                       100732, 69778, 100858,
                                                                       39608, 39708, 71968,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 104008, 0, 3,
                                                                       100858, 69868, 100984,
                                                                       39708, 39808, 72118,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 104218, 0, 3,
                                                                       100984, 69958, 101110,
                                                                       39808, 39908, 72268,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 104428, 0, 3,
                                                                       101110, 70048, 101236,
                                                                       39908, 40008, 72418,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 104638, 0, 3,
                                                                       101236, 70138, 101362,
                                                                       40008, 40108, 72568,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 104848, 0, 3,
                                                                       101488, 70318, 101698,
                                                                       40308, 40458, 72718,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 105163, 0, 3,
                                                                       101698, 70468, 101908,
                                                                       40458, 40608, 72943,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 105478, 0, 3,
                                                                       101908, 70618, 102118,
                                                                       40608, 40758, 73168,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 105793, 0, 3,
                                                                       102118, 70768, 102328,
                                                                       40758, 40908, 73393,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 106108, 0, 3,
                                                                       102328, 70918, 102538,
                                                                       40908, 41058, 73618,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 106423, 0, 3,
                                                                       102538, 71068, 102748,
                                                                       41058, 41208, 73843,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 106738, 0, 3,
                                                                       102748, 71218, 102958,
                                                                       41208, 41358, 74068,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 107053, 0, 3,
                                                                       103168, 71518, 103378,
                                                                       41658, 41808, 74293,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 107368, 0, 3,
                                                                       103378, 71668, 103588,
                                                                       41808, 41958, 74518,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 107683, 0, 3,
                                                                       103588, 71818, 103798,
                                                                       41958, 42108, 74743,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 107998, 0, 3,
                                                                       103798, 71968, 104008,
                                                                       42108, 42258, 74968,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 108313, 0, 3,
                                                                       104008, 72118, 104218,
                                                                       42258, 42408, 75193,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 108628, 0, 3,
                                                                       104218, 72268, 104428,
                                                                       42408, 42558, 75418,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 108943, 0, 3,
                                                                       104428, 72418, 104638,
                                                                       42558, 42708, 75643,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 109258, 0, 3,
                                                                       104848, 72718, 105163,
                                                                       43008, 43218, 75868,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 109699, 0, 3,
                                                                       105163, 72943, 105478,
                                                                       43218, 43428, 76183,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 110140, 0, 3,
                                                                       105478, 73168, 105793,
                                                                       43428, 43638, 76498,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 110581, 0, 3,
                                                                       105793, 73393, 106108,
                                                                       43638, 43848, 76813,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 111022, 0, 3,
                                                                       106108, 73618, 106423,
                                                                       43848, 44058, 77128,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 111463, 0, 3,
                                                                       106423, 73843, 106738,
                                                                       44058, 44268, 77443,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 111904, 0, 3,
                                                                       107053, 74293, 107368,
                                                                       44688, 44898, 77758,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 112345, 0, 3,
                                                                       107368, 74518, 107683,
                                                                       44898, 45108, 78073,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 112786, 0, 3,
                                                                       107683, 74743, 107998,
                                                                       45108, 45318, 78388,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 113227, 0, 3,
                                                                       107998, 74968, 108313,
                                                                       45318, 45528, 78703,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 113668, 0, 3,
                                                                       108313, 75193, 108628,
                                                                       45528, 45738, 79018,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 114109, 0, 3,
                                                                       108628, 75418, 108943,
                                                                       45738, 45948, 79333,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 114550, 0, 3,
                                                                       109258, 75868, 109699,
                                                                       46368, 46648, 79648,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 115138, 0, 3,
                                                                       109699, 76183, 110140,
                                                                       46648, 46928, 80068,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 115726, 0, 3,
                                                                       110140, 76498, 110581,
                                                                       46928, 47208, 80488,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 116314, 0, 3,
                                                                       110581, 76813, 111022,
                                                                       47208, 47488, 80908,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 116902, 0, 3,
                                                                       111022, 77128, 111463,
                                                                       47488, 47768, 81328,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 117490, 0, 3,
                                                                       111904, 77758, 112345,
                                                                       48328, 48608, 81748,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 118078, 0, 3,
                                                                       112345, 78073, 112786,
                                                                       48608, 48888, 82168,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 118666, 0, 3,
                                                                       112786, 78388, 113227,
                                                                       48888, 49168, 82588,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 119254, 0, 3,
                                                                       113227, 78703, 113668,
                                                                       49168, 49448, 83008,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 119842, 0, 3,
                                                                       113668, 79018, 114109,
                                                                       49448, 49728, 83428,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 120430, 0, 3,
                                                                       114550, 79648, 115138,
                                                                       50288, 50648, 83848,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 121186, 0, 3,
                                                                       115138, 80068, 115726,
                                                                       50648, 51008, 84388,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 121942, 0, 3,
                                                                       115726, 80488, 116314,
                                                                       51008, 51368, 84928,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 122698, 0, 3,
                                                                       116314, 80908, 116902,
                                                                       51368, 51728, 85468,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 123454, 0, 3,
                                                                       117490, 81748, 118078,
                                                                       52448, 52808, 86008,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 124210, 0, 3,
                                                                       118078, 82168, 118666,
                                                                       52808, 53168, 86548,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 124966, 0, 3,
                                                                       118666, 82588, 119254,
                                                                       53168, 53528, 87088,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 125722, 0, 3,
                                                                       119254, 83008, 119842,
                                                                       53528, 53888, 87628,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 126478, 0, 3,
                                                                       120430, 83848, 121186,
                                                                       54608, 55058, 88168,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 127423, 0, 3,
                                                                       121186, 84388, 121942,
                                                                       55058, 55508, 88843,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 128368, 0, 3,
                                                                       121942, 84928, 122698,
                                                                       55508, 55958, 89518,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 129313, 0, 3,
                                                                       123454, 86008, 124210,
                                                                       56858, 57308, 90193,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 130258, 0, 3,
                                                                       124210, 86548, 124966,
                                                                       57308, 57758, 90868,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 131203, 0, 3,
                                                                       124966, 87088, 125722,
                                                                       57758, 58208, 91543,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 132148, 0, 3,
                                                                       126478, 88168, 127423,
                                                                       59108, 59658, 92218,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 133303, 0, 3,
                                                                       127423, 88843, 128368,
                                                                       59658, 60208, 93043,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 134458, 0, 3,
                                                                       129313, 90193, 130258,
                                                                       61308, 61858, 93868,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 135613, 0, 3,
                                                                       130258, 90868, 131203,
                                                                       61858, 62408, 94693,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 136768, 0, 3,
                                                                       132148, 92218, 133303,
                                                                       63508, 64168, 95518,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 138154, 0, 3,
                                                                       134458, 93868, 135613,
                                                                       65488, 66148, 96508,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 139540, 114550, 588, ncols);

                    simdfunc::contract_primitives(buffer, 140436, 117490, 588, ncols);

                    simdfunc::contract_primitives(buffer, 141332, 120430, 756, ncols);

                    simdfunc::contract_primitives(buffer, 142484, 123454, 756, ncols);

                    simdfunc::contract_primitives(buffer, 143636, 126478, 945, ncols);

                    simdfunc::contract_primitives(buffer, 145076, 129313, 945, ncols);

                    simdfunc::contract_primitives(buffer, 146516, 132148, 1155, ncols);

                    simdfunc::contract_primitives(buffer, 148276, 134458, 1155, ncols);

                    simdfunc::contract_primitives(buffer, 150036, 136768, 1386, ncols);

                    simdfunc::contract_primitives(buffer, 152148, 138154, 1386, ncols);
                }
            }
        }

        simdtrf::transform_h_inner(buffer, 140128, 139540, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 141024, 140436, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 142088, 141332, 36, 1, nmax);

        simdtrf::transform_h_inner(buffer, 143240, 142484, 36, 1, nmax);

        simdtrf::transform_h_inner(buffer, 144581, 143636, 45, 1, nmax);

        simdtrf::transform_h_inner(buffer, 146021, 145076, 45, 1, nmax);

        simdtrf::transform_h_inner(buffer, 147671, 146516, 55, 1, nmax);

        simdtrf::transform_h_inner(buffer, 149431, 148276, 55, 1, nmax);

        simdtrf::transform_h_inner(buffer, 151422, 150036, 66, 1, nmax);

        simdtrf::transform_h_inner(buffer, 153534, 152148, 66, 1, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 154260, 140128, 142088, 11, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 155184, 141024, 143240, 11, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 156108, 142088, 144581, 11, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 157296, 143240, 146021, 11, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 158484, 144581, 147671, 11, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 159969, 146021, 149431, 11, nmax);

        simdtrf::compute_hrr_pm(buffer, coordinates, 161454, 147671, 151422, 11, nmax);

        simdtrf::compute_hrr_pm(buffer, coordinates, 163269, 149431, 153534, 11, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 165084, 154260, 156108, 11, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 166932, 155184, 157296, 11, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 168780, 156108, 158484, 11, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 171156, 157296, 159969, 11, nmax);

        simdtrf::compute_hrr_dl(buffer, coordinates, 173532, 158484, 161454, 11, nmax);

        simdtrf::compute_hrr_dl(buffer, coordinates, 176502, 159969, 163269, 11, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 179472, 165084, 168780, 11, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 182552, 166932, 171156, 11, nmax);

        simdtrf::compute_hrr_fk(buffer, coordinates, 185632, 168780, 173532, 11, nmax);

        simdtrf::compute_hrr_fk(buffer, coordinates, 189592, 171156, 176502, 11, nmax);

        simdtrf::compute_hrr_gi(buffer, coordinates, 193552, 179472, 185632, 11, nmax);

        simdtrf::compute_hrr_gi(buffer, coordinates, 198172, 182552, 189592, 11, nmax);

        simdtrf::transform_i_inner(buffer, 202792, 198172, 15, 11, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 202792, 143, nmax);

        simdtrf::transform_i_inner(buffer, 202792, 193552, 15, 11, nmax);

        simdtrf::transform_g_outer(values + 1287 * nvalues + n * npairs, nvalues, buffer, 202792,
                                   143, nmax);
    }

    for (size_t m = 0; m < 2574; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
