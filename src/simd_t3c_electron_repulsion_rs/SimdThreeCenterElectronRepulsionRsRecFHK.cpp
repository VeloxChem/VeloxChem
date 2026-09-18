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


#include "SimdThreeCenterElectronRepulsionRsRecFHK.hpp"

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
#include "SimdTransferDH.hpp"
#include "SimdTransferDI.hpp"
#include "SimdTransferFH.hpp"
#include "SimdTransferPH.hpp"
#include "SimdTransferPI.hpp"
#include "SimdTransferPK.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformK.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_fhk_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_fhk_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 209252, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 2310 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 209252, 171572, 12585, dimensions);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2978, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2981, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2984, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2987, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2990, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2993, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2996, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2999, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3002, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3005, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3008, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3011, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3014, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3017, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3020, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3023, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3026, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3029, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3032, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3035, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3038, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3041, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3044, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3047, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3050, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3053, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3056, 3, 34,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3059, 3, 35,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3062, 3, 36,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3065, 3, 37,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3068, 3, 9, 44,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3077, 3, 10, 47,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3086, 3, 11, 50,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3095, 3, 12, 53,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3104, 3, 13, 56,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3113, 3, 14, 59,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3122, 3, 15, 62,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3131, 3, 16, 65,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3140, 3, 17, 68,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3149, 3, 18, 71,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3158, 3, 19, 74,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3167, 3, 20, 77,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3176, 3, 25, 86,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3185, 3, 26, 89,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3194, 3, 27, 92,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3203, 3, 28, 95,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3212, 3, 29, 98,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3221, 3, 30, 101,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3230, 3, 31, 104,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3239, 3, 32, 107,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3248, 3, 33, 110,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3257, 3, 34, 113,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3266, 3, 35, 116,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3275, 3, 36, 119,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3284, 3, 38, 122,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3302, 3, 41, 128,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3320, 3, 44, 134,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3338, 3, 47, 140,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3356, 3, 50, 146,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3374, 3, 53, 152,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3392, 3, 56, 158,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3410, 3, 59, 164,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3428, 3, 62, 170,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3446, 3, 65, 176,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3464, 3, 68, 182,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3482, 3, 71, 188,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3500, 3, 74, 194,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3518, 3, 80, 200,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3536, 3, 83, 206,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3554, 3, 86, 212,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3572, 3, 89, 218,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3590, 3, 92, 224,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3608, 3, 95, 230,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3626, 3, 98, 236,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3644, 3, 101, 242,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3662, 3, 104, 248,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3680, 3, 107, 254,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3698, 3, 110, 260,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3716, 3, 113, 266,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3734, 3, 116, 272,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3752, 3, 122, 278,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3782, 3, 128, 288,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3812, 3, 134, 298,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3842, 3, 140, 308,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3872, 3, 146, 318,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3902, 3, 152, 328,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3932, 3, 158, 338,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3962, 3, 164, 348,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3992, 3, 170, 358,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4022, 3, 176, 368,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4052, 3, 182, 378,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4082, 3, 188, 388,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4112, 3, 200, 398,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4142, 3, 206, 408,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4172, 3, 212, 418,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4202, 3, 218, 428,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4232, 3, 224, 438,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4262, 3, 230, 448,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4292, 3, 236, 458,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4322, 3, 242, 468,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4352, 3, 248, 478,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4382, 3, 254, 488,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4412, 3, 260, 498,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4442, 3, 266, 508,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4472, 3, 278, 518,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4517, 3, 288, 533,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4562, 3, 298, 548,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4607, 3, 308, 563,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4652, 3, 318, 578,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4697, 3, 328, 593,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4742, 3, 338, 608,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4787, 3, 348, 623,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4832, 3, 358, 638,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4877, 3, 368, 653,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4922, 3, 378, 668,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4967, 3, 398, 683,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5012, 3, 408, 698,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5057, 3, 418, 713,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5102, 3, 428, 728,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5147, 3, 438, 743,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5192, 3, 448, 758,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5237, 3, 458, 773,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5282, 3, 468, 788,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5327, 3, 478, 803,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5372, 3, 488, 818,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5417, 3, 498, 833,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5462, 3, 518, 848,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5525, 3, 533, 869,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5588, 3, 548, 890,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5651, 3, 563, 911,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5714, 3, 578, 932,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5777, 3, 593, 953,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5840, 3, 608, 974,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5903, 3, 623, 995,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5966, 3, 638,
                                                                       1016, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6029, 3, 653,
                                                                       1037, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6092, 3, 683,
                                                                       1058, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6155, 3, 698,
                                                                       1079, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6218, 3, 713,
                                                                       1100, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6281, 3, 728,
                                                                       1121, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6344, 3, 743,
                                                                       1142, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6407, 3, 758,
                                                                       1163, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6470, 3, 773,
                                                                       1184, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6533, 3, 788,
                                                                       1205, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6596, 3, 803,
                                                                       1226, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6659, 3, 818,
                                                                       1247, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6722, 3, 848,
                                                                       1268, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6806, 3, 869,
                                                                       1296, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6890, 3, 890,
                                                                       1324, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6974, 3, 911,
                                                                       1352, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7058, 3, 932,
                                                                       1380, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7142, 3, 953,
                                                                       1408, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7226, 3, 974,
                                                                       1436, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7310, 3, 995,
                                                                       1464, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7394, 3, 1016,
                                                                       1492, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7478, 3, 1058,
                                                                       1520, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7562, 3, 1079,
                                                                       1548, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7646, 3, 1100,
                                                                       1576, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7730, 3, 1121,
                                                                       1604, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7814, 3, 1142,
                                                                       1632, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7898, 3, 1163,
                                                                       1660, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7982, 3, 1184,
                                                                       1688, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8066, 3, 1205,
                                                                       1716, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8150, 3, 1226,
                                                                       1744, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8234, 3, 1268,
                                                                       1772, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8342, 3, 1296,
                                                                       1808, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8450, 3, 1324,
                                                                       1844, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8558, 3, 1352,
                                                                       1880, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8666, 3, 1380,
                                                                       1916, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8774, 3, 1408,
                                                                       1952, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8882, 3, 1436,
                                                                       1988, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8990, 3, 1464,
                                                                       2024, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 9098, 3, 1520,
                                                                       2060, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 9206, 3, 1548,
                                                                       2096, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 9314, 3, 1576,
                                                                       2132, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 9422, 3, 1604,
                                                                       2168, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 9530, 3, 1632,
                                                                       2204, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 9638, 3, 1660,
                                                                       2240, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 9746, 3, 1688,
                                                                       2276, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 9854, 3, 1716,
                                                                       2312, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 9962, 3, 1772,
                                                                       2348, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 10097, 3, 1808,
                                                                       2393, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 10232, 3, 1844,
                                                                       2438, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 10367, 3, 1880,
                                                                       2483, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 10502, 3, 1916,
                                                                       2528, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 10637, 3, 1952,
                                                                       2573, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 10772, 3, 1988,
                                                                       2618, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 10907, 3, 2060,
                                                                       2663, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 11042, 3, 2096,
                                                                       2708, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 11177, 3, 2132,
                                                                       2753, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 11312, 3, 2168,
                                                                       2798, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 11447, 3, 2204,
                                                                       2843, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 11582, 3, 2240,
                                                                       2888, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 11717, 3, 2276,
                                                                       2933, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11852, 3, 7, 8,
                                                                       2984, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11858, 3, 8, 9,
                                                                       2987, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11864, 3, 9, 10,
                                                                       2990, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11870, 3, 10, 11,
                                                                       2993, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11876, 3, 11, 12,
                                                                       2996, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11882, 3, 12, 13,
                                                                       2999, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11888, 3, 13, 14,
                                                                       3002, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11894, 3, 14, 15,
                                                                       3005, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11900, 3, 15, 16,
                                                                       3008, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11906, 3, 16, 17,
                                                                       3011, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11912, 3, 17, 18,
                                                                       3014, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11918, 3, 18, 19,
                                                                       3017, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11924, 3, 19, 20,
                                                                       3020, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11930, 3, 23, 24,
                                                                       3029, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11936, 3, 24, 25,
                                                                       3032, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11942, 3, 25, 26,
                                                                       3035, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11948, 3, 26, 27,
                                                                       3038, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11954, 3, 27, 28,
                                                                       3041, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11960, 3, 28, 29,
                                                                       3044, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11966, 3, 29, 30,
                                                                       3047, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11972, 3, 30, 31,
                                                                       3050, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11978, 3, 31, 32,
                                                                       3053, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11984, 3, 32, 33,
                                                                       3056, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11990, 3, 33, 34,
                                                                       3059, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11996, 3, 34, 35,
                                                                       3062, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12002, 3, 35, 36,
                                                                       3065, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12008, 0, 3,
                                                                       11852, 2984, 11858, 3068,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12026, 0, 3,
                                                                       11858, 2987, 11864, 3077,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12044, 0, 3,
                                                                       11864, 2990, 11870, 3086,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12062, 0, 3,
                                                                       11870, 2993, 11876, 3095,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12080, 0, 3,
                                                                       11876, 2996, 11882, 3104,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12098, 0, 3,
                                                                       11882, 2999, 11888, 3113,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12116, 0, 3,
                                                                       11888, 3002, 11894, 3122,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12134, 0, 3,
                                                                       11894, 3005, 11900, 3131,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12152, 0, 3,
                                                                       11900, 3008, 11906, 3140,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12170, 0, 3,
                                                                       11906, 3011, 11912, 3149,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12188, 0, 3,
                                                                       11912, 3014, 11918, 3158,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12206, 0, 3,
                                                                       11918, 3017, 11924, 3167,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12224, 0, 3,
                                                                       11930, 3029, 11936, 3176,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12242, 0, 3,
                                                                       11936, 3032, 11942, 3185,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12260, 0, 3,
                                                                       11942, 3035, 11948, 3194,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12278, 0, 3,
                                                                       11948, 3038, 11954, 3203,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12296, 0, 3,
                                                                       11954, 3041, 11960, 3212,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12314, 0, 3,
                                                                       11960, 3044, 11966, 3221,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12332, 0, 3,
                                                                       11966, 3047, 11972, 3230,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12350, 0, 3,
                                                                       11972, 3050, 11978, 3239,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12368, 0, 3,
                                                                       11978, 3053, 11984, 3248,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12386, 0, 3,
                                                                       11984, 3056, 11990, 3257,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12404, 0, 3,
                                                                       11990, 3059, 11996, 3266,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12422, 0, 3,
                                                                       11996, 3062, 12002, 3275,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 12440, 0, 3,
                                                                       12008, 3068, 12026, 122,
                                                                       128, 3320, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 12476, 0, 3,
                                                                       12026, 3077, 12044, 128,
                                                                       134, 3338, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 12512, 0, 3,
                                                                       12044, 3086, 12062, 134,
                                                                       140, 3356, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 12548, 0, 3,
                                                                       12062, 3095, 12080, 140,
                                                                       146, 3374, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 12584, 0, 3,
                                                                       12080, 3104, 12098, 146,
                                                                       152, 3392, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 12620, 0, 3,
                                                                       12098, 3113, 12116, 152,
                                                                       158, 3410, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 12656, 0, 3,
                                                                       12116, 3122, 12134, 158,
                                                                       164, 3428, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 12692, 0, 3,
                                                                       12134, 3131, 12152, 164,
                                                                       170, 3446, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 12728, 0, 3,
                                                                       12152, 3140, 12170, 170,
                                                                       176, 3464, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 12764, 0, 3,
                                                                       12170, 3149, 12188, 176,
                                                                       182, 3482, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 12800, 0, 3,
                                                                       12188, 3158, 12206, 182,
                                                                       188, 3500, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 12836, 0, 3,
                                                                       12224, 3176, 12242, 200,
                                                                       206, 3554, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 12872, 0, 3,
                                                                       12242, 3185, 12260, 206,
                                                                       212, 3572, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 12908, 0, 3,
                                                                       12260, 3194, 12278, 212,
                                                                       218, 3590, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 12944, 0, 3,
                                                                       12278, 3203, 12296, 218,
                                                                       224, 3608, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 12980, 0, 3,
                                                                       12296, 3212, 12314, 224,
                                                                       230, 3626, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13016, 0, 3,
                                                                       12314, 3221, 12332, 230,
                                                                       236, 3644, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13052, 0, 3,
                                                                       12332, 3230, 12350, 236,
                                                                       242, 3662, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13088, 0, 3,
                                                                       12350, 3239, 12368, 242,
                                                                       248, 3680, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13124, 0, 3,
                                                                       12368, 3248, 12386, 248,
                                                                       254, 3698, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13160, 0, 3,
                                                                       12386, 3257, 12404, 254,
                                                                       260, 3716, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13196, 0, 3,
                                                                       12404, 3266, 12422, 260,
                                                                       266, 3734, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 13232, 0, 3,
                                                                       12440, 3320, 12476, 278,
                                                                       288, 3812, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 13292, 0, 3,
                                                                       12476, 3338, 12512, 288,
                                                                       298, 3842, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 13352, 0, 3,
                                                                       12512, 3356, 12548, 298,
                                                                       308, 3872, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 13412, 0, 3,
                                                                       12548, 3374, 12584, 308,
                                                                       318, 3902, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 13472, 0, 3,
                                                                       12584, 3392, 12620, 318,
                                                                       328, 3932, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 13532, 0, 3,
                                                                       12620, 3410, 12656, 328,
                                                                       338, 3962, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 13592, 0, 3,
                                                                       12656, 3428, 12692, 338,
                                                                       348, 3992, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 13652, 0, 3,
                                                                       12692, 3446, 12728, 348,
                                                                       358, 4022, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 13712, 0, 3,
                                                                       12728, 3464, 12764, 358,
                                                                       368, 4052, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 13772, 0, 3,
                                                                       12764, 3482, 12800, 368,
                                                                       378, 4082, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 13832, 0, 3,
                                                                       12836, 3554, 12872, 398,
                                                                       408, 4172, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 13892, 0, 3,
                                                                       12872, 3572, 12908, 408,
                                                                       418, 4202, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 13952, 0, 3,
                                                                       12908, 3590, 12944, 418,
                                                                       428, 4232, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 14012, 0, 3,
                                                                       12944, 3608, 12980, 428,
                                                                       438, 4262, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 14072, 0, 3,
                                                                       12980, 3626, 13016, 438,
                                                                       448, 4292, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 14132, 0, 3,
                                                                       13016, 3644, 13052, 448,
                                                                       458, 4322, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 14192, 0, 3,
                                                                       13052, 3662, 13088, 458,
                                                                       468, 4352, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 14252, 0, 3,
                                                                       13088, 3680, 13124, 468,
                                                                       478, 4382, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 14312, 0, 3,
                                                                       13124, 3698, 13160, 478,
                                                                       488, 4412, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 14372, 0, 3,
                                                                       13160, 3716, 13196, 488,
                                                                       498, 4442, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 14432, 0, 3,
                                                                       13232, 3812, 13292, 518,
                                                                       533, 4562, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 14522, 0, 3,
                                                                       13292, 3842, 13352, 533,
                                                                       548, 4607, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 14612, 0, 3,
                                                                       13352, 3872, 13412, 548,
                                                                       563, 4652, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 14702, 0, 3,
                                                                       13412, 3902, 13472, 563,
                                                                       578, 4697, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 14792, 0, 3,
                                                                       13472, 3932, 13532, 578,
                                                                       593, 4742, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 14882, 0, 3,
                                                                       13532, 3962, 13592, 593,
                                                                       608, 4787, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 14972, 0, 3,
                                                                       13592, 3992, 13652, 608,
                                                                       623, 4832, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 15062, 0, 3,
                                                                       13652, 4022, 13712, 623,
                                                                       638, 4877, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 15152, 0, 3,
                                                                       13712, 4052, 13772, 638,
                                                                       653, 4922, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 15242, 0, 3,
                                                                       13832, 4172, 13892, 683,
                                                                       698, 5057, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 15332, 0, 3,
                                                                       13892, 4202, 13952, 698,
                                                                       713, 5102, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 15422, 0, 3,
                                                                       13952, 4232, 14012, 713,
                                                                       728, 5147, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 15512, 0, 3,
                                                                       14012, 4262, 14072, 728,
                                                                       743, 5192, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 15602, 0, 3,
                                                                       14072, 4292, 14132, 743,
                                                                       758, 5237, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 15692, 0, 3,
                                                                       14132, 4322, 14192, 758,
                                                                       773, 5282, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 15782, 0, 3,
                                                                       14192, 4352, 14252, 773,
                                                                       788, 5327, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 15872, 0, 3,
                                                                       14252, 4382, 14312, 788,
                                                                       803, 5372, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 15962, 0, 3,
                                                                       14312, 4412, 14372, 803,
                                                                       818, 5417, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 16052, 0, 3,
                                                                       14432, 4562, 14522, 848,
                                                                       869, 5588, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 16178, 0, 3,
                                                                       14522, 4607, 14612, 869,
                                                                       890, 5651, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 16304, 0, 3,
                                                                       14612, 4652, 14702, 890,
                                                                       911, 5714, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 16430, 0, 3,
                                                                       14702, 4697, 14792, 911,
                                                                       932, 5777, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 16556, 0, 3,
                                                                       14792, 4742, 14882, 932,
                                                                       953, 5840, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 16682, 0, 3,
                                                                       14882, 4787, 14972, 953,
                                                                       974, 5903, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 16808, 0, 3,
                                                                       14972, 4832, 15062, 974,
                                                                       995, 5966, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 16934, 0, 3,
                                                                       15062, 4877, 15152, 995,
                                                                       1016, 6029, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 17060, 0, 3,
                                                                       15242, 5057, 15332, 1058,
                                                                       1079, 6218, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 17186, 0, 3,
                                                                       15332, 5102, 15422, 1079,
                                                                       1100, 6281, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 17312, 0, 3,
                                                                       15422, 5147, 15512, 1100,
                                                                       1121, 6344, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 17438, 0, 3,
                                                                       15512, 5192, 15602, 1121,
                                                                       1142, 6407, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 17564, 0, 3,
                                                                       15602, 5237, 15692, 1142,
                                                                       1163, 6470, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 17690, 0, 3,
                                                                       15692, 5282, 15782, 1163,
                                                                       1184, 6533, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 17816, 0, 3,
                                                                       15782, 5327, 15872, 1184,
                                                                       1205, 6596, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 17942, 0, 3,
                                                                       15872, 5372, 15962, 1205,
                                                                       1226, 6659, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 18068, 0, 3,
                                                                       16052, 5588, 16178, 1268,
                                                                       1296, 6890, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 18236, 0, 3,
                                                                       16178, 5651, 16304, 1296,
                                                                       1324, 6974, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 18404, 0, 3,
                                                                       16304, 5714, 16430, 1324,
                                                                       1352, 7058, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 18572, 0, 3,
                                                                       16430, 5777, 16556, 1352,
                                                                       1380, 7142, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 18740, 0, 3,
                                                                       16556, 5840, 16682, 1380,
                                                                       1408, 7226, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 18908, 0, 3,
                                                                       16682, 5903, 16808, 1408,
                                                                       1436, 7310, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 19076, 0, 3,
                                                                       16808, 5966, 16934, 1436,
                                                                       1464, 7394, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 19244, 0, 3,
                                                                       17060, 6218, 17186, 1520,
                                                                       1548, 7646, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 19412, 0, 3,
                                                                       17186, 6281, 17312, 1548,
                                                                       1576, 7730, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 19580, 0, 3,
                                                                       17312, 6344, 17438, 1576,
                                                                       1604, 7814, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 19748, 0, 3,
                                                                       17438, 6407, 17564, 1604,
                                                                       1632, 7898, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 19916, 0, 3,
                                                                       17564, 6470, 17690, 1632,
                                                                       1660, 7982, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 20084, 0, 3,
                                                                       17690, 6533, 17816, 1660,
                                                                       1688, 8066, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 20252, 0, 3,
                                                                       17816, 6596, 17942, 1688,
                                                                       1716, 8150, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 20420, 0, 3,
                                                                       18068, 6890, 18236, 1772,
                                                                       1808, 8450, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 20636, 0, 3,
                                                                       18236, 6974, 18404, 1808,
                                                                       1844, 8558, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 20852, 0, 3,
                                                                       18404, 7058, 18572, 1844,
                                                                       1880, 8666, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 21068, 0, 3,
                                                                       18572, 7142, 18740, 1880,
                                                                       1916, 8774, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 21284, 0, 3,
                                                                       18740, 7226, 18908, 1916,
                                                                       1952, 8882, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 21500, 0, 3,
                                                                       18908, 7310, 19076, 1952,
                                                                       1988, 8990, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 21716, 0, 3,
                                                                       19244, 7646, 19412, 2060,
                                                                       2096, 9314, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 21932, 0, 3,
                                                                       19412, 7730, 19580, 2096,
                                                                       2132, 9422, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 22148, 0, 3,
                                                                       19580, 7814, 19748, 2132,
                                                                       2168, 9530, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 22364, 0, 3,
                                                                       19748, 7898, 19916, 2168,
                                                                       2204, 9638, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 22580, 0, 3,
                                                                       19916, 7982, 20084, 2204,
                                                                       2240, 9746, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 22796, 0, 3,
                                                                       20084, 8066, 20252, 2240,
                                                                       2276, 9854, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 23012, 0, 3,
                                                                       20420, 8450, 20636, 2348,
                                                                       2393, 10232, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 23282, 0, 3,
                                                                       20636, 8558, 20852, 2393,
                                                                       2438, 10367, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 23552, 0, 3,
                                                                       20852, 8666, 21068, 2438,
                                                                       2483, 10502, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 23822, 0, 3,
                                                                       21068, 8774, 21284, 2483,
                                                                       2528, 10637, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 24092, 0, 3,
                                                                       21284, 8882, 21500, 2528,
                                                                       2573, 10772, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 24362, 0, 3,
                                                                       21716, 9314, 21932, 2663,
                                                                       2708, 11177, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 24632, 0, 3,
                                                                       21932, 9422, 22148, 2708,
                                                                       2753, 11312, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 24902, 0, 3,
                                                                       22148, 9530, 22364, 2753,
                                                                       2798, 11447, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 25172, 0, 3,
                                                                       22364, 9638, 22580, 2798,
                                                                       2843, 11582, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 25442, 0, 3,
                                                                       22580, 9746, 22796, 2843,
                                                                       2888, 11717, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25712, 3, 2978,
                                                                       2981, 11852, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25722, 3, 2981,
                                                                       2984, 11858, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25732, 3, 2984,
                                                                       2987, 11864, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25742, 3, 2987,
                                                                       2990, 11870, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25752, 3, 2990,
                                                                       2993, 11876, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25762, 3, 2993,
                                                                       2996, 11882, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25772, 3, 2996,
                                                                       2999, 11888, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25782, 3, 2999,
                                                                       3002, 11894, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25792, 3, 3002,
                                                                       3005, 11900, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25802, 3, 3005,
                                                                       3008, 11906, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25812, 3, 3008,
                                                                       3011, 11912, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25822, 3, 3011,
                                                                       3014, 11918, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25832, 3, 3014,
                                                                       3017, 11924, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25842, 3, 3023,
                                                                       3026, 11930, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25852, 3, 3026,
                                                                       3029, 11936, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25862, 3, 3029,
                                                                       3032, 11942, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25872, 3, 3032,
                                                                       3035, 11948, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25882, 3, 3035,
                                                                       3038, 11954, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25892, 3, 3038,
                                                                       3041, 11960, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25902, 3, 3041,
                                                                       3044, 11966, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25912, 3, 3044,
                                                                       3047, 11972, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25922, 3, 3047,
                                                                       3050, 11978, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25932, 3, 3050,
                                                                       3053, 11984, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25942, 3, 3053,
                                                                       3056, 11990, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25952, 3, 3056,
                                                                       3059, 11996, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25962, 3, 3059,
                                                                       3062, 12002, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 25972, 0, 3,
                                                                       25712, 11852, 25722,
                                                                       12008, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 26002, 0, 3,
                                                                       25722, 11858, 25732,
                                                                       12026, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 26032, 0, 3,
                                                                       25732, 11864, 25742,
                                                                       12044, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 26062, 0, 3,
                                                                       25742, 11870, 25752,
                                                                       12062, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 26092, 0, 3,
                                                                       25752, 11876, 25762,
                                                                       12080, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 26122, 0, 3,
                                                                       25762, 11882, 25772,
                                                                       12098, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 26152, 0, 3,
                                                                       25772, 11888, 25782,
                                                                       12116, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 26182, 0, 3,
                                                                       25782, 11894, 25792,
                                                                       12134, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 26212, 0, 3,
                                                                       25792, 11900, 25802,
                                                                       12152, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 26242, 0, 3,
                                                                       25802, 11906, 25812,
                                                                       12170, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 26272, 0, 3,
                                                                       25812, 11912, 25822,
                                                                       12188, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 26302, 0, 3,
                                                                       25822, 11918, 25832,
                                                                       12206, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 26332, 0, 3,
                                                                       25842, 11930, 25852,
                                                                       12224, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 26362, 0, 3,
                                                                       25852, 11936, 25862,
                                                                       12242, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 26392, 0, 3,
                                                                       25862, 11942, 25872,
                                                                       12260, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 26422, 0, 3,
                                                                       25872, 11948, 25882,
                                                                       12278, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 26452, 0, 3,
                                                                       25882, 11954, 25892,
                                                                       12296, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 26482, 0, 3,
                                                                       25892, 11960, 25902,
                                                                       12314, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 26512, 0, 3,
                                                                       25902, 11966, 25912,
                                                                       12332, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 26542, 0, 3,
                                                                       25912, 11972, 25922,
                                                                       12350, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 26572, 0, 3,
                                                                       25922, 11978, 25932,
                                                                       12368, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 26602, 0, 3,
                                                                       25932, 11984, 25942,
                                                                       12386, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 26632, 0, 3,
                                                                       25942, 11990, 25952,
                                                                       12404, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 26662, 0, 3,
                                                                       25952, 11996, 25962,
                                                                       12422, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 26692, 0, 3,
                                                                       25972, 12008, 26002, 3284,
                                                                       3302, 12440, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 26752, 0, 3,
                                                                       26002, 12026, 26032, 3302,
                                                                       3320, 12476, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 26812, 0, 3,
                                                                       26032, 12044, 26062, 3320,
                                                                       3338, 12512, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 26872, 0, 3,
                                                                       26062, 12062, 26092, 3338,
                                                                       3356, 12548, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 26932, 0, 3,
                                                                       26092, 12080, 26122, 3356,
                                                                       3374, 12584, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 26992, 0, 3,
                                                                       26122, 12098, 26152, 3374,
                                                                       3392, 12620, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 27052, 0, 3,
                                                                       26152, 12116, 26182, 3392,
                                                                       3410, 12656, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 27112, 0, 3,
                                                                       26182, 12134, 26212, 3410,
                                                                       3428, 12692, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 27172, 0, 3,
                                                                       26212, 12152, 26242, 3428,
                                                                       3446, 12728, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 27232, 0, 3,
                                                                       26242, 12170, 26272, 3446,
                                                                       3464, 12764, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 27292, 0, 3,
                                                                       26272, 12188, 26302, 3464,
                                                                       3482, 12800, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 27352, 0, 3,
                                                                       26332, 12224, 26362, 3518,
                                                                       3536, 12836, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 27412, 0, 3,
                                                                       26362, 12242, 26392, 3536,
                                                                       3554, 12872, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 27472, 0, 3,
                                                                       26392, 12260, 26422, 3554,
                                                                       3572, 12908, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 27532, 0, 3,
                                                                       26422, 12278, 26452, 3572,
                                                                       3590, 12944, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 27592, 0, 3,
                                                                       26452, 12296, 26482, 3590,
                                                                       3608, 12980, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 27652, 0, 3,
                                                                       26482, 12314, 26512, 3608,
                                                                       3626, 13016, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 27712, 0, 3,
                                                                       26512, 12332, 26542, 3626,
                                                                       3644, 13052, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 27772, 0, 3,
                                                                       26542, 12350, 26572, 3644,
                                                                       3662, 13088, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 27832, 0, 3,
                                                                       26572, 12368, 26602, 3662,
                                                                       3680, 13124, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 27892, 0, 3,
                                                                       26602, 12386, 26632, 3680,
                                                                       3698, 13160, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 27952, 0, 3,
                                                                       26632, 12404, 26662, 3698,
                                                                       3716, 13196, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 28012, 0, 3,
                                                                       26692, 12440, 26752, 3752,
                                                                       3782, 13232, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 28112, 0, 3,
                                                                       26752, 12476, 26812, 3782,
                                                                       3812, 13292, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 28212, 0, 3,
                                                                       26812, 12512, 26872, 3812,
                                                                       3842, 13352, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 28312, 0, 3,
                                                                       26872, 12548, 26932, 3842,
                                                                       3872, 13412, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 28412, 0, 3,
                                                                       26932, 12584, 26992, 3872,
                                                                       3902, 13472, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 28512, 0, 3,
                                                                       26992, 12620, 27052, 3902,
                                                                       3932, 13532, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 28612, 0, 3,
                                                                       27052, 12656, 27112, 3932,
                                                                       3962, 13592, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 28712, 0, 3,
                                                                       27112, 12692, 27172, 3962,
                                                                       3992, 13652, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 28812, 0, 3,
                                                                       27172, 12728, 27232, 3992,
                                                                       4022, 13712, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 28912, 0, 3,
                                                                       27232, 12764, 27292, 4022,
                                                                       4052, 13772, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 29012, 0, 3,
                                                                       27352, 12836, 27412, 4112,
                                                                       4142, 13832, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 29112, 0, 3,
                                                                       27412, 12872, 27472, 4142,
                                                                       4172, 13892, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 29212, 0, 3,
                                                                       27472, 12908, 27532, 4172,
                                                                       4202, 13952, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 29312, 0, 3,
                                                                       27532, 12944, 27592, 4202,
                                                                       4232, 14012, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 29412, 0, 3,
                                                                       27592, 12980, 27652, 4232,
                                                                       4262, 14072, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 29512, 0, 3,
                                                                       27652, 13016, 27712, 4262,
                                                                       4292, 14132, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 29612, 0, 3,
                                                                       27712, 13052, 27772, 4292,
                                                                       4322, 14192, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 29712, 0, 3,
                                                                       27772, 13088, 27832, 4322,
                                                                       4352, 14252, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 29812, 0, 3,
                                                                       27832, 13124, 27892, 4352,
                                                                       4382, 14312, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 29912, 0, 3,
                                                                       27892, 13160, 27952, 4382,
                                                                       4412, 14372, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 30012, 0, 3,
                                                                       28012, 13232, 28112, 4472,
                                                                       4517, 14432, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 30162, 0, 3,
                                                                       28112, 13292, 28212, 4517,
                                                                       4562, 14522, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 30312, 0, 3,
                                                                       28212, 13352, 28312, 4562,
                                                                       4607, 14612, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 30462, 0, 3,
                                                                       28312, 13412, 28412, 4607,
                                                                       4652, 14702, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 30612, 0, 3,
                                                                       28412, 13472, 28512, 4652,
                                                                       4697, 14792, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 30762, 0, 3,
                                                                       28512, 13532, 28612, 4697,
                                                                       4742, 14882, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 30912, 0, 3,
                                                                       28612, 13592, 28712, 4742,
                                                                       4787, 14972, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 31062, 0, 3,
                                                                       28712, 13652, 28812, 4787,
                                                                       4832, 15062, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 31212, 0, 3,
                                                                       28812, 13712, 28912, 4832,
                                                                       4877, 15152, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 31362, 0, 3,
                                                                       29012, 13832, 29112, 4967,
                                                                       5012, 15242, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 31512, 0, 3,
                                                                       29112, 13892, 29212, 5012,
                                                                       5057, 15332, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 31662, 0, 3,
                                                                       29212, 13952, 29312, 5057,
                                                                       5102, 15422, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 31812, 0, 3,
                                                                       29312, 14012, 29412, 5102,
                                                                       5147, 15512, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 31962, 0, 3,
                                                                       29412, 14072, 29512, 5147,
                                                                       5192, 15602, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 32112, 0, 3,
                                                                       29512, 14132, 29612, 5192,
                                                                       5237, 15692, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 32262, 0, 3,
                                                                       29612, 14192, 29712, 5237,
                                                                       5282, 15782, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 32412, 0, 3,
                                                                       29712, 14252, 29812, 5282,
                                                                       5327, 15872, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 32562, 0, 3,
                                                                       29812, 14312, 29912, 5327,
                                                                       5372, 15962, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 32712, 0, 3,
                                                                       30012, 14432, 30162, 5462,
                                                                       5525, 16052, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 32922, 0, 3,
                                                                       30162, 14522, 30312, 5525,
                                                                       5588, 16178, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 33132, 0, 3,
                                                                       30312, 14612, 30462, 5588,
                                                                       5651, 16304, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 33342, 0, 3,
                                                                       30462, 14702, 30612, 5651,
                                                                       5714, 16430, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 33552, 0, 3,
                                                                       30612, 14792, 30762, 5714,
                                                                       5777, 16556, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 33762, 0, 3,
                                                                       30762, 14882, 30912, 5777,
                                                                       5840, 16682, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 33972, 0, 3,
                                                                       30912, 14972, 31062, 5840,
                                                                       5903, 16808, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 34182, 0, 3,
                                                                       31062, 15062, 31212, 5903,
                                                                       5966, 16934, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 34392, 0, 3,
                                                                       31362, 15242, 31512, 6092,
                                                                       6155, 17060, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 34602, 0, 3,
                                                                       31512, 15332, 31662, 6155,
                                                                       6218, 17186, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 34812, 0, 3,
                                                                       31662, 15422, 31812, 6218,
                                                                       6281, 17312, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 35022, 0, 3,
                                                                       31812, 15512, 31962, 6281,
                                                                       6344, 17438, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 35232, 0, 3,
                                                                       31962, 15602, 32112, 6344,
                                                                       6407, 17564, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 35442, 0, 3,
                                                                       32112, 15692, 32262, 6407,
                                                                       6470, 17690, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 35652, 0, 3,
                                                                       32262, 15782, 32412, 6470,
                                                                       6533, 17816, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 35862, 0, 3,
                                                                       32412, 15872, 32562, 6533,
                                                                       6596, 17942, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 36072, 0, 3,
                                                                       32712, 16052, 32922, 6722,
                                                                       6806, 18068, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 36352, 0, 3,
                                                                       32922, 16178, 33132, 6806,
                                                                       6890, 18236, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 36632, 0, 3,
                                                                       33132, 16304, 33342, 6890,
                                                                       6974, 18404, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 36912, 0, 3,
                                                                       33342, 16430, 33552, 6974,
                                                                       7058, 18572, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 37192, 0, 3,
                                                                       33552, 16556, 33762, 7058,
                                                                       7142, 18740, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 37472, 0, 3,
                                                                       33762, 16682, 33972, 7142,
                                                                       7226, 18908, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 37752, 0, 3,
                                                                       33972, 16808, 34182, 7226,
                                                                       7310, 19076, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 38032, 0, 3,
                                                                       34392, 17060, 34602, 7478,
                                                                       7562, 19244, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 38312, 0, 3,
                                                                       34602, 17186, 34812, 7562,
                                                                       7646, 19412, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 38592, 0, 3,
                                                                       34812, 17312, 35022, 7646,
                                                                       7730, 19580, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 38872, 0, 3,
                                                                       35022, 17438, 35232, 7730,
                                                                       7814, 19748, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 39152, 0, 3,
                                                                       35232, 17564, 35442, 7814,
                                                                       7898, 19916, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 39432, 0, 3,
                                                                       35442, 17690, 35652, 7898,
                                                                       7982, 20084, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 39712, 0, 3,
                                                                       35652, 17816, 35862, 7982,
                                                                       8066, 20252, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 39992, 0, 3,
                                                                       36072, 18068, 36352, 8234,
                                                                       8342, 20420, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 40352, 0, 3,
                                                                       36352, 18236, 36632, 8342,
                                                                       8450, 20636, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 40712, 0, 3,
                                                                       36632, 18404, 36912, 8450,
                                                                       8558, 20852, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 41072, 0, 3,
                                                                       36912, 18572, 37192, 8558,
                                                                       8666, 21068, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 41432, 0, 3,
                                                                       37192, 18740, 37472, 8666,
                                                                       8774, 21284, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 41792, 0, 3,
                                                                       37472, 18908, 37752, 8774,
                                                                       8882, 21500, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 42152, 0, 3,
                                                                       38032, 19244, 38312, 9098,
                                                                       9206, 21716, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 42512, 0, 3,
                                                                       38312, 19412, 38592, 9206,
                                                                       9314, 21932, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 42872, 0, 3,
                                                                       38592, 19580, 38872, 9314,
                                                                       9422, 22148, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 43232, 0, 3,
                                                                       38872, 19748, 39152, 9422,
                                                                       9530, 22364, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 43592, 0, 3,
                                                                       39152, 19916, 39432, 9530,
                                                                       9638, 22580, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 43952, 0, 3,
                                                                       39432, 20084, 39712, 9638,
                                                                       9746, 22796, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 44312, 0, 3,
                                                                       39992, 20420, 40352, 9962,
                                                                       10097, 23012, ncols,
                                                                       gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 44762, 0, 3,
                                                                       40352, 20636, 40712,
                                                                       10097, 10232, 23282,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 45212, 0, 3,
                                                                       40712, 20852, 41072,
                                                                       10232, 10367, 23552,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 45662, 0, 3,
                                                                       41072, 21068, 41432,
                                                                       10367, 10502, 23822,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 46112, 0, 3,
                                                                       41432, 21284, 41792,
                                                                       10502, 10637, 24092,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 46562, 0, 3,
                                                                       42152, 21716, 42512,
                                                                       10907, 11042, 24362,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 47012, 0, 3,
                                                                       42512, 21932, 42872,
                                                                       11042, 11177, 24632,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 47462, 0, 3,
                                                                       42872, 22148, 43232,
                                                                       11177, 11312, 24902,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 47912, 0, 3,
                                                                       43232, 22364, 43592,
                                                                       11312, 11447, 25172,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 48362, 0, 3,
                                                                       43592, 22580, 43952,
                                                                       11447, 11582, 25442,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 48812, 3, 11852,
                                                                       11858, 25732, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 48827, 3, 11858,
                                                                       11864, 25742, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 48842, 3, 11864,
                                                                       11870, 25752, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 48857, 3, 11870,
                                                                       11876, 25762, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 48872, 3, 11876,
                                                                       11882, 25772, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 48887, 3, 11882,
                                                                       11888, 25782, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 48902, 3, 11888,
                                                                       11894, 25792, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 48917, 3, 11894,
                                                                       11900, 25802, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 48932, 3, 11900,
                                                                       11906, 25812, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 48947, 3, 11906,
                                                                       11912, 25822, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 48962, 3, 11912,
                                                                       11918, 25832, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 48977, 3, 11930,
                                                                       11936, 25862, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 48992, 3, 11936,
                                                                       11942, 25872, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49007, 3, 11942,
                                                                       11948, 25882, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49022, 3, 11948,
                                                                       11954, 25892, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49037, 3, 11954,
                                                                       11960, 25902, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49052, 3, 11960,
                                                                       11966, 25912, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49067, 3, 11966,
                                                                       11972, 25922, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49082, 3, 11972,
                                                                       11978, 25932, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49097, 3, 11978,
                                                                       11984, 25942, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49112, 3, 11984,
                                                                       11990, 25952, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49127, 3, 11990,
                                                                       11996, 25962, ncols,
                                                                       gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 49142, 0, 3,
                                                                       48812, 25732, 48827,
                                                                       12008, 12026, 26032,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 49187, 0, 3,
                                                                       48827, 25742, 48842,
                                                                       12026, 12044, 26062,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 49232, 0, 3,
                                                                       48842, 25752, 48857,
                                                                       12044, 12062, 26092,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 49277, 0, 3,
                                                                       48857, 25762, 48872,
                                                                       12062, 12080, 26122,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 49322, 0, 3,
                                                                       48872, 25772, 48887,
                                                                       12080, 12098, 26152,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 49367, 0, 3,
                                                                       48887, 25782, 48902,
                                                                       12098, 12116, 26182,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 49412, 0, 3,
                                                                       48902, 25792, 48917,
                                                                       12116, 12134, 26212,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 49457, 0, 3,
                                                                       48917, 25802, 48932,
                                                                       12134, 12152, 26242,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 49502, 0, 3,
                                                                       48932, 25812, 48947,
                                                                       12152, 12170, 26272,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 49547, 0, 3,
                                                                       48947, 25822, 48962,
                                                                       12170, 12188, 26302,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 49592, 0, 3,
                                                                       48977, 25862, 48992,
                                                                       12224, 12242, 26392,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 49637, 0, 3,
                                                                       48992, 25872, 49007,
                                                                       12242, 12260, 26422,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 49682, 0, 3,
                                                                       49007, 25882, 49022,
                                                                       12260, 12278, 26452,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 49727, 0, 3,
                                                                       49022, 25892, 49037,
                                                                       12278, 12296, 26482,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 49772, 0, 3,
                                                                       49037, 25902, 49052,
                                                                       12296, 12314, 26512,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 49817, 0, 3,
                                                                       49052, 25912, 49067,
                                                                       12314, 12332, 26542,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 49862, 0, 3,
                                                                       49067, 25922, 49082,
                                                                       12332, 12350, 26572,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 49907, 0, 3,
                                                                       49082, 25932, 49097,
                                                                       12350, 12368, 26602,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 49952, 0, 3,
                                                                       49097, 25942, 49112,
                                                                       12368, 12386, 26632,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 49997, 0, 3,
                                                                       49112, 25952, 49127,
                                                                       12386, 12404, 26662,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 50042, 0, 3,
                                                                       49142, 26032, 49187,
                                                                       12440, 12476, 26812,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 50132, 0, 3,
                                                                       49187, 26062, 49232,
                                                                       12476, 12512, 26872,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 50222, 0, 3,
                                                                       49232, 26092, 49277,
                                                                       12512, 12548, 26932,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 50312, 0, 3,
                                                                       49277, 26122, 49322,
                                                                       12548, 12584, 26992,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 50402, 0, 3,
                                                                       49322, 26152, 49367,
                                                                       12584, 12620, 27052,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 50492, 0, 3,
                                                                       49367, 26182, 49412,
                                                                       12620, 12656, 27112,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 50582, 0, 3,
                                                                       49412, 26212, 49457,
                                                                       12656, 12692, 27172,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 50672, 0, 3,
                                                                       49457, 26242, 49502,
                                                                       12692, 12728, 27232,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 50762, 0, 3,
                                                                       49502, 26272, 49547,
                                                                       12728, 12764, 27292,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 50852, 0, 3,
                                                                       49592, 26392, 49637,
                                                                       12836, 12872, 27472,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 50942, 0, 3,
                                                                       49637, 26422, 49682,
                                                                       12872, 12908, 27532,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 51032, 0, 3,
                                                                       49682, 26452, 49727,
                                                                       12908, 12944, 27592,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 51122, 0, 3,
                                                                       49727, 26482, 49772,
                                                                       12944, 12980, 27652,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 51212, 0, 3,
                                                                       49772, 26512, 49817,
                                                                       12980, 13016, 27712,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 51302, 0, 3,
                                                                       49817, 26542, 49862,
                                                                       13016, 13052, 27772,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 51392, 0, 3,
                                                                       49862, 26572, 49907,
                                                                       13052, 13088, 27832,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 51482, 0, 3,
                                                                       49907, 26602, 49952,
                                                                       13088, 13124, 27892,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 51572, 0, 3,
                                                                       49952, 26632, 49997,
                                                                       13124, 13160, 27952,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 51662, 0, 3,
                                                                       50042, 26812, 50132,
                                                                       13232, 13292, 28212,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 51812, 0, 3,
                                                                       50132, 26872, 50222,
                                                                       13292, 13352, 28312,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 51962, 0, 3,
                                                                       50222, 26932, 50312,
                                                                       13352, 13412, 28412,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 52112, 0, 3,
                                                                       50312, 26992, 50402,
                                                                       13412, 13472, 28512,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 52262, 0, 3,
                                                                       50402, 27052, 50492,
                                                                       13472, 13532, 28612,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 52412, 0, 3,
                                                                       50492, 27112, 50582,
                                                                       13532, 13592, 28712,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 52562, 0, 3,
                                                                       50582, 27172, 50672,
                                                                       13592, 13652, 28812,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 52712, 0, 3,
                                                                       50672, 27232, 50762,
                                                                       13652, 13712, 28912,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 52862, 0, 3,
                                                                       50852, 27472, 50942,
                                                                       13832, 13892, 29212,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 53012, 0, 3,
                                                                       50942, 27532, 51032,
                                                                       13892, 13952, 29312,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 53162, 0, 3,
                                                                       51032, 27592, 51122,
                                                                       13952, 14012, 29412,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 53312, 0, 3,
                                                                       51122, 27652, 51212,
                                                                       14012, 14072, 29512,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 53462, 0, 3,
                                                                       51212, 27712, 51302,
                                                                       14072, 14132, 29612,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 53612, 0, 3,
                                                                       51302, 27772, 51392,
                                                                       14132, 14192, 29712,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 53762, 0, 3,
                                                                       51392, 27832, 51482,
                                                                       14192, 14252, 29812,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 53912, 0, 3,
                                                                       51482, 27892, 51572,
                                                                       14252, 14312, 29912,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 54062, 0, 3,
                                                                       51662, 28212, 51812,
                                                                       14432, 14522, 30312,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 54287, 0, 3,
                                                                       51812, 28312, 51962,
                                                                       14522, 14612, 30462,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 54512, 0, 3,
                                                                       51962, 28412, 52112,
                                                                       14612, 14702, 30612,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 54737, 0, 3,
                                                                       52112, 28512, 52262,
                                                                       14702, 14792, 30762,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 54962, 0, 3,
                                                                       52262, 28612, 52412,
                                                                       14792, 14882, 30912,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 55187, 0, 3,
                                                                       52412, 28712, 52562,
                                                                       14882, 14972, 31062,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 55412, 0, 3,
                                                                       52562, 28812, 52712,
                                                                       14972, 15062, 31212,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 55637, 0, 3,
                                                                       52862, 29212, 53012,
                                                                       15242, 15332, 31662,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 55862, 0, 3,
                                                                       53012, 29312, 53162,
                                                                       15332, 15422, 31812,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 56087, 0, 3,
                                                                       53162, 29412, 53312,
                                                                       15422, 15512, 31962,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 56312, 0, 3,
                                                                       53312, 29512, 53462,
                                                                       15512, 15602, 32112,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 56537, 0, 3,
                                                                       53462, 29612, 53612,
                                                                       15602, 15692, 32262,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 56762, 0, 3,
                                                                       53612, 29712, 53762,
                                                                       15692, 15782, 32412,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 56987, 0, 3,
                                                                       53762, 29812, 53912,
                                                                       15782, 15872, 32562,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 57212, 0, 3,
                                                                       54062, 30312, 54287,
                                                                       16052, 16178, 33132,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 57527, 0, 3,
                                                                       54287, 30462, 54512,
                                                                       16178, 16304, 33342,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 57842, 0, 3,
                                                                       54512, 30612, 54737,
                                                                       16304, 16430, 33552,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 58157, 0, 3,
                                                                       54737, 30762, 54962,
                                                                       16430, 16556, 33762,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 58472, 0, 3,
                                                                       54962, 30912, 55187,
                                                                       16556, 16682, 33972,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 58787, 0, 3,
                                                                       55187, 31062, 55412,
                                                                       16682, 16808, 34182,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 59102, 0, 3,
                                                                       55637, 31662, 55862,
                                                                       17060, 17186, 34812,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 59417, 0, 3,
                                                                       55862, 31812, 56087,
                                                                       17186, 17312, 35022,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 59732, 0, 3,
                                                                       56087, 31962, 56312,
                                                                       17312, 17438, 35232,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 60047, 0, 3,
                                                                       56312, 32112, 56537,
                                                                       17438, 17564, 35442,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 60362, 0, 3,
                                                                       56537, 32262, 56762,
                                                                       17564, 17690, 35652,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 60677, 0, 3,
                                                                       56762, 32412, 56987,
                                                                       17690, 17816, 35862,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 60992, 0, 3,
                                                                       57212, 33132, 57527,
                                                                       18068, 18236, 36632,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 61412, 0, 3,
                                                                       57527, 33342, 57842,
                                                                       18236, 18404, 36912,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 61832, 0, 3,
                                                                       57842, 33552, 58157,
                                                                       18404, 18572, 37192,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 62252, 0, 3,
                                                                       58157, 33762, 58472,
                                                                       18572, 18740, 37472,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 62672, 0, 3,
                                                                       58472, 33972, 58787,
                                                                       18740, 18908, 37752,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 63092, 0, 3,
                                                                       59102, 34812, 59417,
                                                                       19244, 19412, 38592,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 63512, 0, 3,
                                                                       59417, 35022, 59732,
                                                                       19412, 19580, 38872,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 63932, 0, 3,
                                                                       59732, 35232, 60047,
                                                                       19580, 19748, 39152,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 64352, 0, 3,
                                                                       60047, 35442, 60362,
                                                                       19748, 19916, 39432,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 64772, 0, 3,
                                                                       60362, 35652, 60677,
                                                                       19916, 20084, 39712,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 65192, 0, 3,
                                                                       60992, 36632, 61412,
                                                                       20420, 20636, 40712,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 65732, 0, 3,
                                                                       61412, 36912, 61832,
                                                                       20636, 20852, 41072,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 66272, 0, 3,
                                                                       61832, 37192, 62252,
                                                                       20852, 21068, 41432,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 66812, 0, 3,
                                                                       62252, 37472, 62672,
                                                                       21068, 21284, 41792,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 67352, 0, 3,
                                                                       63092, 38592, 63512,
                                                                       21716, 21932, 42872,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 67892, 0, 3,
                                                                       63512, 38872, 63932,
                                                                       21932, 22148, 43232,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 68432, 0, 3,
                                                                       63932, 39152, 64352,
                                                                       22148, 22364, 43592,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 68972, 0, 3,
                                                                       64352, 39432, 64772,
                                                                       22364, 22580, 43952,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 69512, 0, 3,
                                                                       65192, 40712, 65732,
                                                                       23012, 23282, 45212,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 70187, 0, 3,
                                                                       65732, 41072, 66272,
                                                                       23282, 23552, 45662,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 70862, 0, 3,
                                                                       66272, 41432, 66812,
                                                                       23552, 23822, 46112,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 71537, 0, 3,
                                                                       67352, 42872, 67892,
                                                                       24362, 24632, 47462,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 72212, 0, 3,
                                                                       67892, 43232, 68432,
                                                                       24632, 24902, 47912,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 72887, 0, 3,
                                                                       68432, 43592, 68972,
                                                                       24902, 25172, 48362,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 73562, 3, 25712,
                                                                       25722, 48812, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 73583, 3, 25722,
                                                                       25732, 48827, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 73604, 3, 25732,
                                                                       25742, 48842, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 73625, 3, 25742,
                                                                       25752, 48857, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 73646, 3, 25752,
                                                                       25762, 48872, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 73667, 3, 25762,
                                                                       25772, 48887, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 73688, 3, 25772,
                                                                       25782, 48902, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 73709, 3, 25782,
                                                                       25792, 48917, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 73730, 3, 25792,
                                                                       25802, 48932, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 73751, 3, 25802,
                                                                       25812, 48947, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 73772, 3, 25812,
                                                                       25822, 48962, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 73793, 3, 25842,
                                                                       25852, 48977, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 73814, 3, 25852,
                                                                       25862, 48992, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 73835, 3, 25862,
                                                                       25872, 49007, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 73856, 3, 25872,
                                                                       25882, 49022, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 73877, 3, 25882,
                                                                       25892, 49037, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 73898, 3, 25892,
                                                                       25902, 49052, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 73919, 3, 25902,
                                                                       25912, 49067, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 73940, 3, 25912,
                                                                       25922, 49082, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 73961, 3, 25922,
                                                                       25932, 49097, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 73982, 3, 25932,
                                                                       25942, 49112, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 74003, 3, 25942,
                                                                       25952, 49127, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 74024, 0, 3,
                                                                       73562, 48812, 73583,
                                                                       25972, 26002, 49142,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 74087, 0, 3,
                                                                       73583, 48827, 73604,
                                                                       26002, 26032, 49187,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 74150, 0, 3,
                                                                       73604, 48842, 73625,
                                                                       26032, 26062, 49232,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 74213, 0, 3,
                                                                       73625, 48857, 73646,
                                                                       26062, 26092, 49277,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 74276, 0, 3,
                                                                       73646, 48872, 73667,
                                                                       26092, 26122, 49322,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 74339, 0, 3,
                                                                       73667, 48887, 73688,
                                                                       26122, 26152, 49367,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 74402, 0, 3,
                                                                       73688, 48902, 73709,
                                                                       26152, 26182, 49412,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 74465, 0, 3,
                                                                       73709, 48917, 73730,
                                                                       26182, 26212, 49457,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 74528, 0, 3,
                                                                       73730, 48932, 73751,
                                                                       26212, 26242, 49502,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 74591, 0, 3,
                                                                       73751, 48947, 73772,
                                                                       26242, 26272, 49547,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 74654, 0, 3,
                                                                       73793, 48977, 73814,
                                                                       26332, 26362, 49592,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 74717, 0, 3,
                                                                       73814, 48992, 73835,
                                                                       26362, 26392, 49637,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 74780, 0, 3,
                                                                       73835, 49007, 73856,
                                                                       26392, 26422, 49682,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 74843, 0, 3,
                                                                       73856, 49022, 73877,
                                                                       26422, 26452, 49727,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 74906, 0, 3,
                                                                       73877, 49037, 73898,
                                                                       26452, 26482, 49772,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 74969, 0, 3,
                                                                       73898, 49052, 73919,
                                                                       26482, 26512, 49817,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 75032, 0, 3,
                                                                       73919, 49067, 73940,
                                                                       26512, 26542, 49862,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 75095, 0, 3,
                                                                       73940, 49082, 73961,
                                                                       26542, 26572, 49907,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 75158, 0, 3,
                                                                       73961, 49097, 73982,
                                                                       26572, 26602, 49952,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 75221, 0, 3,
                                                                       73982, 49112, 74003,
                                                                       26602, 26632, 49997,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 75284, 0, 3,
                                                                       74024, 49142, 74087,
                                                                       26692, 26752, 50042,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 75410, 0, 3,
                                                                       74087, 49187, 74150,
                                                                       26752, 26812, 50132,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 75536, 0, 3,
                                                                       74150, 49232, 74213,
                                                                       26812, 26872, 50222,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 75662, 0, 3,
                                                                       74213, 49277, 74276,
                                                                       26872, 26932, 50312,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 75788, 0, 3,
                                                                       74276, 49322, 74339,
                                                                       26932, 26992, 50402,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 75914, 0, 3,
                                                                       74339, 49367, 74402,
                                                                       26992, 27052, 50492,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 76040, 0, 3,
                                                                       74402, 49412, 74465,
                                                                       27052, 27112, 50582,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 76166, 0, 3,
                                                                       74465, 49457, 74528,
                                                                       27112, 27172, 50672,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 76292, 0, 3,
                                                                       74528, 49502, 74591,
                                                                       27172, 27232, 50762,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 76418, 0, 3,
                                                                       74654, 49592, 74717,
                                                                       27352, 27412, 50852,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 76544, 0, 3,
                                                                       74717, 49637, 74780,
                                                                       27412, 27472, 50942,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 76670, 0, 3,
                                                                       74780, 49682, 74843,
                                                                       27472, 27532, 51032,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 76796, 0, 3,
                                                                       74843, 49727, 74906,
                                                                       27532, 27592, 51122,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 76922, 0, 3,
                                                                       74906, 49772, 74969,
                                                                       27592, 27652, 51212,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 77048, 0, 3,
                                                                       74969, 49817, 75032,
                                                                       27652, 27712, 51302,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 77174, 0, 3,
                                                                       75032, 49862, 75095,
                                                                       27712, 27772, 51392,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 77300, 0, 3,
                                                                       75095, 49907, 75158,
                                                                       27772, 27832, 51482,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 77426, 0, 3,
                                                                       75158, 49952, 75221,
                                                                       27832, 27892, 51572,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 77552, 0, 3,
                                                                       75284, 50042, 75410,
                                                                       28012, 28112, 51662,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 77762, 0, 3,
                                                                       75410, 50132, 75536,
                                                                       28112, 28212, 51812,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 77972, 0, 3,
                                                                       75536, 50222, 75662,
                                                                       28212, 28312, 51962,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 78182, 0, 3,
                                                                       75662, 50312, 75788,
                                                                       28312, 28412, 52112,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 78392, 0, 3,
                                                                       75788, 50402, 75914,
                                                                       28412, 28512, 52262,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 78602, 0, 3,
                                                                       75914, 50492, 76040,
                                                                       28512, 28612, 52412,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 78812, 0, 3,
                                                                       76040, 50582, 76166,
                                                                       28612, 28712, 52562,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 79022, 0, 3,
                                                                       76166, 50672, 76292,
                                                                       28712, 28812, 52712,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 79232, 0, 3,
                                                                       76418, 50852, 76544,
                                                                       29012, 29112, 52862,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 79442, 0, 3,
                                                                       76544, 50942, 76670,
                                                                       29112, 29212, 53012,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 79652, 0, 3,
                                                                       76670, 51032, 76796,
                                                                       29212, 29312, 53162,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 79862, 0, 3,
                                                                       76796, 51122, 76922,
                                                                       29312, 29412, 53312,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 80072, 0, 3,
                                                                       76922, 51212, 77048,
                                                                       29412, 29512, 53462,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 80282, 0, 3,
                                                                       77048, 51302, 77174,
                                                                       29512, 29612, 53612,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 80492, 0, 3,
                                                                       77174, 51392, 77300,
                                                                       29612, 29712, 53762,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 80702, 0, 3,
                                                                       77300, 51482, 77426,
                                                                       29712, 29812, 53912,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 80912, 0, 3,
                                                                       77552, 51662, 77762,
                                                                       30012, 30162, 54062,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 81227, 0, 3,
                                                                       77762, 51812, 77972,
                                                                       30162, 30312, 54287,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 81542, 0, 3,
                                                                       77972, 51962, 78182,
                                                                       30312, 30462, 54512,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 81857, 0, 3,
                                                                       78182, 52112, 78392,
                                                                       30462, 30612, 54737,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 82172, 0, 3,
                                                                       78392, 52262, 78602,
                                                                       30612, 30762, 54962,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 82487, 0, 3,
                                                                       78602, 52412, 78812,
                                                                       30762, 30912, 55187,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 82802, 0, 3,
                                                                       78812, 52562, 79022,
                                                                       30912, 31062, 55412,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 83117, 0, 3,
                                                                       79232, 52862, 79442,
                                                                       31362, 31512, 55637,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 83432, 0, 3,
                                                                       79442, 53012, 79652,
                                                                       31512, 31662, 55862,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 83747, 0, 3,
                                                                       79652, 53162, 79862,
                                                                       31662, 31812, 56087,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 84062, 0, 3,
                                                                       79862, 53312, 80072,
                                                                       31812, 31962, 56312,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 84377, 0, 3,
                                                                       80072, 53462, 80282,
                                                                       31962, 32112, 56537,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 84692, 0, 3,
                                                                       80282, 53612, 80492,
                                                                       32112, 32262, 56762,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 85007, 0, 3,
                                                                       80492, 53762, 80702,
                                                                       32262, 32412, 56987,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 85322, 0, 3,
                                                                       80912, 54062, 81227,
                                                                       32712, 32922, 57212,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 85763, 0, 3,
                                                                       81227, 54287, 81542,
                                                                       32922, 33132, 57527,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 86204, 0, 3,
                                                                       81542, 54512, 81857,
                                                                       33132, 33342, 57842,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 86645, 0, 3,
                                                                       81857, 54737, 82172,
                                                                       33342, 33552, 58157,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 87086, 0, 3,
                                                                       82172, 54962, 82487,
                                                                       33552, 33762, 58472,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 87527, 0, 3,
                                                                       82487, 55187, 82802,
                                                                       33762, 33972, 58787,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 87968, 0, 3,
                                                                       83117, 55637, 83432,
                                                                       34392, 34602, 59102,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 88409, 0, 3,
                                                                       83432, 55862, 83747,
                                                                       34602, 34812, 59417,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 88850, 0, 3,
                                                                       83747, 56087, 84062,
                                                                       34812, 35022, 59732,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 89291, 0, 3,
                                                                       84062, 56312, 84377,
                                                                       35022, 35232, 60047,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 89732, 0, 3,
                                                                       84377, 56537, 84692,
                                                                       35232, 35442, 60362,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 90173, 0, 3,
                                                                       84692, 56762, 85007,
                                                                       35442, 35652, 60677,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 90614, 0, 3,
                                                                       85322, 57212, 85763,
                                                                       36072, 36352, 60992,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 91202, 0, 3,
                                                                       85763, 57527, 86204,
                                                                       36352, 36632, 61412,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 91790, 0, 3,
                                                                       86204, 57842, 86645,
                                                                       36632, 36912, 61832,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 92378, 0, 3,
                                                                       86645, 58157, 87086,
                                                                       36912, 37192, 62252,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 92966, 0, 3,
                                                                       87086, 58472, 87527,
                                                                       37192, 37472, 62672,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 93554, 0, 3,
                                                                       87968, 59102, 88409,
                                                                       38032, 38312, 63092,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 94142, 0, 3,
                                                                       88409, 59417, 88850,
                                                                       38312, 38592, 63512,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 94730, 0, 3,
                                                                       88850, 59732, 89291,
                                                                       38592, 38872, 63932,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 95318, 0, 3,
                                                                       89291, 60047, 89732,
                                                                       38872, 39152, 64352,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 95906, 0, 3,
                                                                       89732, 60362, 90173,
                                                                       39152, 39432, 64772,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 96494, 0, 3,
                                                                       90614, 60992, 91202,
                                                                       39992, 40352, 65192,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 97250, 0, 3,
                                                                       91202, 61412, 91790,
                                                                       40352, 40712, 65732,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 98006, 0, 3,
                                                                       91790, 61832, 92378,
                                                                       40712, 41072, 66272,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 98762, 0, 3,
                                                                       92378, 62252, 92966,
                                                                       41072, 41432, 66812,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 99518, 0, 3,
                                                                       93554, 63092, 94142,
                                                                       42152, 42512, 67352,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 100274, 0, 3,
                                                                       94142, 63512, 94730,
                                                                       42512, 42872, 67892,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 101030, 0, 3,
                                                                       94730, 63932, 95318,
                                                                       42872, 43232, 68432,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 101786, 0, 3,
                                                                       95318, 64352, 95906,
                                                                       43232, 43592, 68972,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 102542, 0, 3,
                                                                       96494, 65192, 97250,
                                                                       44312, 44762, 69512,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 103487, 0, 3,
                                                                       97250, 65732, 98006,
                                                                       44762, 45212, 70187,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 104432, 0, 3,
                                                                       98006, 66272, 98762,
                                                                       45212, 45662, 70862,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 105377, 0, 3,
                                                                       99518, 67352, 100274,
                                                                       46562, 47012, 71537,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 106322, 0, 3,
                                                                       100274, 67892, 101030,
                                                                       47012, 47462, 72212,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 107267, 0, 3,
                                                                       101030, 68432, 101786,
                                                                       47462, 47912, 72887,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 108212, 3, 48812,
                                                                       48827, 73604, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 108240, 3, 48827,
                                                                       48842, 73625, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 108268, 3, 48842,
                                                                       48857, 73646, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 108296, 3, 48857,
                                                                       48872, 73667, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 108324, 3, 48872,
                                                                       48887, 73688, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 108352, 3, 48887,
                                                                       48902, 73709, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 108380, 3, 48902,
                                                                       48917, 73730, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 108408, 3, 48917,
                                                                       48932, 73751, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 108436, 3, 48932,
                                                                       48947, 73772, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 108464, 3, 48977,
                                                                       48992, 73835, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 108492, 3, 48992,
                                                                       49007, 73856, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 108520, 3, 49007,
                                                                       49022, 73877, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 108548, 3, 49022,
                                                                       49037, 73898, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 108576, 3, 49037,
                                                                       49052, 73919, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 108604, 3, 49052,
                                                                       49067, 73940, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 108632, 3, 49067,
                                                                       49082, 73961, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 108660, 3, 49082,
                                                                       49097, 73982, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 108688, 3, 49097,
                                                                       49112, 74003, ncols,
                                                                       gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 108716, 0, 3,
                                                                       108212, 73604, 108240,
                                                                       49142, 49187, 74150,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 108800, 0, 3,
                                                                       108240, 73625, 108268,
                                                                       49187, 49232, 74213,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 108884, 0, 3,
                                                                       108268, 73646, 108296,
                                                                       49232, 49277, 74276,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 108968, 0, 3,
                                                                       108296, 73667, 108324,
                                                                       49277, 49322, 74339,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 109052, 0, 3,
                                                                       108324, 73688, 108352,
                                                                       49322, 49367, 74402,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 109136, 0, 3,
                                                                       108352, 73709, 108380,
                                                                       49367, 49412, 74465,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 109220, 0, 3,
                                                                       108380, 73730, 108408,
                                                                       49412, 49457, 74528,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 109304, 0, 3,
                                                                       108408, 73751, 108436,
                                                                       49457, 49502, 74591,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 109388, 0, 3,
                                                                       108464, 73835, 108492,
                                                                       49592, 49637, 74780,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 109472, 0, 3,
                                                                       108492, 73856, 108520,
                                                                       49637, 49682, 74843,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 109556, 0, 3,
                                                                       108520, 73877, 108548,
                                                                       49682, 49727, 74906,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 109640, 0, 3,
                                                                       108548, 73898, 108576,
                                                                       49727, 49772, 74969,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 109724, 0, 3,
                                                                       108576, 73919, 108604,
                                                                       49772, 49817, 75032,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 109808, 0, 3,
                                                                       108604, 73940, 108632,
                                                                       49817, 49862, 75095,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 109892, 0, 3,
                                                                       108632, 73961, 108660,
                                                                       49862, 49907, 75158,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 109976, 0, 3,
                                                                       108660, 73982, 108688,
                                                                       49907, 49952, 75221,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 110060, 0, 3,
                                                                       108716, 74150, 108800,
                                                                       50042, 50132, 75536,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 110228, 0, 3,
                                                                       108800, 74213, 108884,
                                                                       50132, 50222, 75662,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 110396, 0, 3,
                                                                       108884, 74276, 108968,
                                                                       50222, 50312, 75788,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 110564, 0, 3,
                                                                       108968, 74339, 109052,
                                                                       50312, 50402, 75914,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 110732, 0, 3,
                                                                       109052, 74402, 109136,
                                                                       50402, 50492, 76040,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 110900, 0, 3,
                                                                       109136, 74465, 109220,
                                                                       50492, 50582, 76166,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 111068, 0, 3,
                                                                       109220, 74528, 109304,
                                                                       50582, 50672, 76292,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 111236, 0, 3,
                                                                       109388, 74780, 109472,
                                                                       50852, 50942, 76670,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 111404, 0, 3,
                                                                       109472, 74843, 109556,
                                                                       50942, 51032, 76796,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 111572, 0, 3,
                                                                       109556, 74906, 109640,
                                                                       51032, 51122, 76922,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 111740, 0, 3,
                                                                       109640, 74969, 109724,
                                                                       51122, 51212, 77048,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 111908, 0, 3,
                                                                       109724, 75032, 109808,
                                                                       51212, 51302, 77174,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 112076, 0, 3,
                                                                       109808, 75095, 109892,
                                                                       51302, 51392, 77300,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 112244, 0, 3,
                                                                       109892, 75158, 109976,
                                                                       51392, 51482, 77426,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 112412, 0, 3,
                                                                       110060, 75536, 110228,
                                                                       51662, 51812, 77972,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 112692, 0, 3,
                                                                       110228, 75662, 110396,
                                                                       51812, 51962, 78182,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 112972, 0, 3,
                                                                       110396, 75788, 110564,
                                                                       51962, 52112, 78392,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 113252, 0, 3,
                                                                       110564, 75914, 110732,
                                                                       52112, 52262, 78602,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 113532, 0, 3,
                                                                       110732, 76040, 110900,
                                                                       52262, 52412, 78812,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 113812, 0, 3,
                                                                       110900, 76166, 111068,
                                                                       52412, 52562, 79022,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 114092, 0, 3,
                                                                       111236, 76670, 111404,
                                                                       52862, 53012, 79652,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 114372, 0, 3,
                                                                       111404, 76796, 111572,
                                                                       53012, 53162, 79862,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 114652, 0, 3,
                                                                       111572, 76922, 111740,
                                                                       53162, 53312, 80072,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 114932, 0, 3,
                                                                       111740, 77048, 111908,
                                                                       53312, 53462, 80282,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 115212, 0, 3,
                                                                       111908, 77174, 112076,
                                                                       53462, 53612, 80492,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 115492, 0, 3,
                                                                       112076, 77300, 112244,
                                                                       53612, 53762, 80702,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 115772, 0, 3,
                                                                       112412, 77972, 112692,
                                                                       54062, 54287, 81542,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 116192, 0, 3,
                                                                       112692, 78182, 112972,
                                                                       54287, 54512, 81857,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 116612, 0, 3,
                                                                       112972, 78392, 113252,
                                                                       54512, 54737, 82172,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 117032, 0, 3,
                                                                       113252, 78602, 113532,
                                                                       54737, 54962, 82487,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 117452, 0, 3,
                                                                       113532, 78812, 113812,
                                                                       54962, 55187, 82802,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 117872, 0, 3,
                                                                       114092, 79652, 114372,
                                                                       55637, 55862, 83747,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 118292, 0, 3,
                                                                       114372, 79862, 114652,
                                                                       55862, 56087, 84062,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 118712, 0, 3,
                                                                       114652, 80072, 114932,
                                                                       56087, 56312, 84377,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 119132, 0, 3,
                                                                       114932, 80282, 115212,
                                                                       56312, 56537, 84692,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 119552, 0, 3,
                                                                       115212, 80492, 115492,
                                                                       56537, 56762, 85007,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 119972, 0, 3,
                                                                       115772, 81542, 116192,
                                                                       57212, 57527, 86204,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 120560, 0, 3,
                                                                       116192, 81857, 116612,
                                                                       57527, 57842, 86645,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 121148, 0, 3,
                                                                       116612, 82172, 117032,
                                                                       57842, 58157, 87086,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 121736, 0, 3,
                                                                       117032, 82487, 117452,
                                                                       58157, 58472, 87527,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 122324, 0, 3,
                                                                       117872, 83747, 118292,
                                                                       59102, 59417, 88850,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 122912, 0, 3,
                                                                       118292, 84062, 118712,
                                                                       59417, 59732, 89291,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 123500, 0, 3,
                                                                       118712, 84377, 119132,
                                                                       59732, 60047, 89732,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 124088, 0, 3,
                                                                       119132, 84692, 119552,
                                                                       60047, 60362, 90173,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 124676, 0, 3,
                                                                       119972, 86204, 120560,
                                                                       60992, 61412, 91790,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 125460, 0, 3,
                                                                       120560, 86645, 121148,
                                                                       61412, 61832, 92378,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 126244, 0, 3,
                                                                       121148, 87086, 121736,
                                                                       61832, 62252, 92966,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 127028, 0, 3,
                                                                       122324, 88850, 122912,
                                                                       63092, 63512, 94730,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 127812, 0, 3,
                                                                       122912, 89291, 123500,
                                                                       63512, 63932, 95318,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 128596, 0, 3,
                                                                       123500, 89732, 124088,
                                                                       63932, 64352, 95906,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 129380, 0, 3,
                                                                       124676, 91790, 125460,
                                                                       65192, 65732, 98006,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 130388, 0, 3,
                                                                       125460, 92378, 126244,
                                                                       65732, 66272, 98762,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 131396, 0, 3,
                                                                       127028, 94730, 127812,
                                                                       67352, 67892, 101030,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 132404, 0, 3,
                                                                       127812, 95318, 128596,
                                                                       67892, 68432, 101786,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 133412, 0, 3,
                                                                       129380, 98006, 130388,
                                                                       69512, 70187, 104432,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 134672, 0, 3,
                                                                       131396, 101030, 132404,
                                                                       71537, 72212, 107267,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 135932, 3, 73562,
                                                                       73583, 108212, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 135968, 3, 73583,
                                                                       73604, 108240, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 136004, 3, 73604,
                                                                       73625, 108268, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 136040, 3, 73625,
                                                                       73646, 108296, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 136076, 3, 73646,
                                                                       73667, 108324, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 136112, 3, 73667,
                                                                       73688, 108352, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 136148, 3, 73688,
                                                                       73709, 108380, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 136184, 3, 73709,
                                                                       73730, 108408, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 136220, 3, 73730,
                                                                       73751, 108436, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 136256, 3, 73793,
                                                                       73814, 108464, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 136292, 3, 73814,
                                                                       73835, 108492, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 136328, 3, 73835,
                                                                       73856, 108520, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 136364, 3, 73856,
                                                                       73877, 108548, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 136400, 3, 73877,
                                                                       73898, 108576, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 136436, 3, 73898,
                                                                       73919, 108604, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 136472, 3, 73919,
                                                                       73940, 108632, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 136508, 3, 73940,
                                                                       73961, 108660, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 136544, 3, 73961,
                                                                       73982, 108688, ncols,
                                                                       gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 136580, 0, 3,
                                                                       135932, 108212, 135968,
                                                                       74024, 74087, 108716,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 136688, 0, 3,
                                                                       135968, 108240, 136004,
                                                                       74087, 74150, 108800,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 136796, 0, 3,
                                                                       136004, 108268, 136040,
                                                                       74150, 74213, 108884,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 136904, 0, 3,
                                                                       136040, 108296, 136076,
                                                                       74213, 74276, 108968,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 137012, 0, 3,
                                                                       136076, 108324, 136112,
                                                                       74276, 74339, 109052,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 137120, 0, 3,
                                                                       136112, 108352, 136148,
                                                                       74339, 74402, 109136,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 137228, 0, 3,
                                                                       136148, 108380, 136184,
                                                                       74402, 74465, 109220,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 137336, 0, 3,
                                                                       136184, 108408, 136220,
                                                                       74465, 74528, 109304,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 137444, 0, 3,
                                                                       136256, 108464, 136292,
                                                                       74654, 74717, 109388,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 137552, 0, 3,
                                                                       136292, 108492, 136328,
                                                                       74717, 74780, 109472,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 137660, 0, 3,
                                                                       136328, 108520, 136364,
                                                                       74780, 74843, 109556,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 137768, 0, 3,
                                                                       136364, 108548, 136400,
                                                                       74843, 74906, 109640,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 137876, 0, 3,
                                                                       136400, 108576, 136436,
                                                                       74906, 74969, 109724,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 137984, 0, 3,
                                                                       136436, 108604, 136472,
                                                                       74969, 75032, 109808,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 138092, 0, 3,
                                                                       136472, 108632, 136508,
                                                                       75032, 75095, 109892,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 138200, 0, 3,
                                                                       136508, 108660, 136544,
                                                                       75095, 75158, 109976,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 138308, 0, 3,
                                                                       136580, 108716, 136688,
                                                                       75284, 75410, 110060,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 138524, 0, 3,
                                                                       136688, 108800, 136796,
                                                                       75410, 75536, 110228,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 138740, 0, 3,
                                                                       136796, 108884, 136904,
                                                                       75536, 75662, 110396,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 138956, 0, 3,
                                                                       136904, 108968, 137012,
                                                                       75662, 75788, 110564,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 139172, 0, 3,
                                                                       137012, 109052, 137120,
                                                                       75788, 75914, 110732,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 139388, 0, 3,
                                                                       137120, 109136, 137228,
                                                                       75914, 76040, 110900,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 139604, 0, 3,
                                                                       137228, 109220, 137336,
                                                                       76040, 76166, 111068,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 139820, 0, 3,
                                                                       137444, 109388, 137552,
                                                                       76418, 76544, 111236,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 140036, 0, 3,
                                                                       137552, 109472, 137660,
                                                                       76544, 76670, 111404,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 140252, 0, 3,
                                                                       137660, 109556, 137768,
                                                                       76670, 76796, 111572,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 140468, 0, 3,
                                                                       137768, 109640, 137876,
                                                                       76796, 76922, 111740,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 140684, 0, 3,
                                                                       137876, 109724, 137984,
                                                                       76922, 77048, 111908,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 140900, 0, 3,
                                                                       137984, 109808, 138092,
                                                                       77048, 77174, 112076,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 141116, 0, 3,
                                                                       138092, 109892, 138200,
                                                                       77174, 77300, 112244,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 141332, 0, 3,
                                                                       138308, 110060, 138524,
                                                                       77552, 77762, 112412,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 141692, 0, 3,
                                                                       138524, 110228, 138740,
                                                                       77762, 77972, 112692,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 142052, 0, 3,
                                                                       138740, 110396, 138956,
                                                                       77972, 78182, 112972,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 142412, 0, 3,
                                                                       138956, 110564, 139172,
                                                                       78182, 78392, 113252,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 142772, 0, 3,
                                                                       139172, 110732, 139388,
                                                                       78392, 78602, 113532,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 143132, 0, 3,
                                                                       139388, 110900, 139604,
                                                                       78602, 78812, 113812,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 143492, 0, 3,
                                                                       139820, 111236, 140036,
                                                                       79232, 79442, 114092,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 143852, 0, 3,
                                                                       140036, 111404, 140252,
                                                                       79442, 79652, 114372,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 144212, 0, 3,
                                                                       140252, 111572, 140468,
                                                                       79652, 79862, 114652,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 144572, 0, 3,
                                                                       140468, 111740, 140684,
                                                                       79862, 80072, 114932,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 144932, 0, 3,
                                                                       140684, 111908, 140900,
                                                                       80072, 80282, 115212,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 145292, 0, 3,
                                                                       140900, 112076, 141116,
                                                                       80282, 80492, 115492,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 145652, 0, 3,
                                                                       141332, 112412, 141692,
                                                                       80912, 81227, 115772,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 146192, 0, 3,
                                                                       141692, 112692, 142052,
                                                                       81227, 81542, 116192,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 146732, 0, 3,
                                                                       142052, 112972, 142412,
                                                                       81542, 81857, 116612,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 147272, 0, 3,
                                                                       142412, 113252, 142772,
                                                                       81857, 82172, 117032,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 147812, 0, 3,
                                                                       142772, 113532, 143132,
                                                                       82172, 82487, 117452,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 148352, 0, 3,
                                                                       143492, 114092, 143852,
                                                                       83117, 83432, 117872,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 148892, 0, 3,
                                                                       143852, 114372, 144212,
                                                                       83432, 83747, 118292,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 149432, 0, 3,
                                                                       144212, 114652, 144572,
                                                                       83747, 84062, 118712,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 149972, 0, 3,
                                                                       144572, 114932, 144932,
                                                                       84062, 84377, 119132,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 150512, 0, 3,
                                                                       144932, 115212, 145292,
                                                                       84377, 84692, 119552,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 151052, 0, 3,
                                                                       145652, 115772, 146192,
                                                                       85322, 85763, 119972,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 151808, 0, 3,
                                                                       146192, 116192, 146732,
                                                                       85763, 86204, 120560,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 152564, 0, 3,
                                                                       146732, 116612, 147272,
                                                                       86204, 86645, 121148,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 153320, 0, 3,
                                                                       147272, 117032, 147812,
                                                                       86645, 87086, 121736,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 154076, 0, 3,
                                                                       148352, 117872, 148892,
                                                                       87968, 88409, 122324,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 154832, 0, 3,
                                                                       148892, 118292, 149432,
                                                                       88409, 88850, 122912,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 155588, 0, 3,
                                                                       149432, 118712, 149972,
                                                                       88850, 89291, 123500,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 156344, 0, 3,
                                                                       149972, 119132, 150512,
                                                                       89291, 89732, 124088,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 157100, 0, 3,
                                                                       151052, 119972, 151808,
                                                                       90614, 91202, 124676,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 158108, 0, 3,
                                                                       151808, 120560, 152564,
                                                                       91202, 91790, 125460,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 159116, 0, 3,
                                                                       152564, 121148, 153320,
                                                                       91790, 92378, 126244,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 160124, 0, 3,
                                                                       154076, 122324, 154832,
                                                                       93554, 94142, 127028,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 161132, 0, 3,
                                                                       154832, 122912, 155588,
                                                                       94142, 94730, 127812,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 162140, 0, 3,
                                                                       155588, 123500, 156344,
                                                                       94730, 95318, 128596,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 163148, 0, 3,
                                                                       157100, 124676, 158108,
                                                                       96494, 97250, 129380,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 164444, 0, 3,
                                                                       158108, 125460, 159116,
                                                                       97250, 98006, 130388,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 165740, 0, 3,
                                                                       160124, 127028, 161132,
                                                                       99518, 100274, 131396,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 167036, 0, 3,
                                                                       161132, 127812, 162140,
                                                                       100274, 101030, 132404,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 168332, 0, 3,
                                                                       163148, 129380, 164444,
                                                                       102542, 103487, 133412,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 169952, 0, 3,
                                                                       165740, 131396, 167036,
                                                                       105377, 106322, 134672,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 171572, 151052, 756, ncols);

                    simdfunc::contract_primitives(buffer, 172643, 154076, 756, ncols);

                    simdfunc::contract_primitives(buffer, 173714, 157100, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 175142, 160124, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 176570, 163148, 1296, ncols);

                    simdfunc::contract_primitives(buffer, 178406, 165740, 1296, ncols);

                    simdfunc::contract_primitives(buffer, 180242, 168332, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 182537, 169952, 1620, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 172328, 171572, 21, 1, nmax);

        simdtrf::transform_k_inner(buffer, 173399, 172643, 21, 1, nmax);

        simdtrf::transform_k_inner(buffer, 174722, 173714, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 176150, 175142, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 177866, 176570, 36, 1, nmax);

        simdtrf::transform_k_inner(buffer, 179702, 178406, 36, 1, nmax);

        simdtrf::transform_k_inner(buffer, 181862, 180242, 45, 1, nmax);

        simdtrf::transform_k_inner(buffer, 184157, 182537, 45, 1, nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 184832, 172328, 174722, 15, nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 185777, 173399, 176150, 15, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 186722, 174722, 177866, 15, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 187982, 176150, 179702, 15, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 189242, 177866, 181862, 15, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 190862, 179702, 184157, 15, nmax);

        simdtrf::compute_hrr_dh(buffer, coordinates, 192482, 184832, 186722, 15, nmax);

        simdtrf::compute_hrr_dh(buffer, coordinates, 194372, 185777, 187982, 15, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 196262, 186722, 189242, 15, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 198782, 187982, 190862, 15, nmax);

        simdtrf::compute_hrr_fh(buffer, coordinates, 201302, 192482, 196262, 15, nmax);

        simdtrf::compute_hrr_fh(buffer, coordinates, 204452, 194372, 198782, 15, nmax);

        simdtrf::transform_h_inner(buffer, 207602, 204452, 10, 15, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 207602, 165, nmax);

        simdtrf::transform_h_inner(buffer, 207602, 201302, 10, 15, nmax);

        simdtrf::transform_f_outer(values + 1155 * nvalues + n * npairs, nvalues, buffer, 207602,
                                   165, nmax);
    }

    for (size_t m = 0; m < 2310; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
