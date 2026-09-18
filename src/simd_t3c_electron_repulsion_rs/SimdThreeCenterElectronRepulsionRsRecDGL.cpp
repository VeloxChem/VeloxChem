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


#include "SimdThreeCenterElectronRepulsionRsRecDGL.hpp"

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
#include "SimdTransferDG.hpp"
#include "SimdTransferPG.hpp"
#include "SimdTransferPH.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformL.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_dgl_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_dgl_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 128574, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1530 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 128574, 112988, 7460, dimensions);

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

                    simdfunc::compute_full_t3c_erf_boys_function(buffer, coordinates, 6, 3, 14,
                                                                 ncols, fj, i * nprim_b + j, fq,
                                                                 omega);

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 22, 3, 14,
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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1772, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1775, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1778, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1781, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1784, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1787, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1790, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1793, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1796, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1799, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1802, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1805, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1808, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1811, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1814, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1817, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1820, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1823, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1826, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1829, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1832, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1835, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1838, 3, 34,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1841, 3, 35,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1844, 3, 36,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1847, 3, 37,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1850, 3, 9, 44,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1859, 3, 10, 47,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1868, 3, 11, 50,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1877, 3, 12, 53,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1886, 3, 13, 56,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1895, 3, 14, 59,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1904, 3, 15, 62,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1913, 3, 16, 65,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1922, 3, 17, 68,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1931, 3, 18, 71,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1940, 3, 19, 74,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1949, 3, 20, 77,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1958, 3, 25, 86,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1967, 3, 26, 89,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1976, 3, 27, 92,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1985, 3, 28, 95,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1994, 3, 29, 98,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2003, 3, 30, 101,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2012, 3, 31, 104,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2021, 3, 32, 107,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2030, 3, 33, 110,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2039, 3, 34, 113,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2048, 3, 35, 116,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2057, 3, 36, 119,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2066, 3, 44, 134,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2084, 3, 47, 140,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2102, 3, 50, 146,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2120, 3, 53, 152,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2138, 3, 56, 158,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2156, 3, 59, 164,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2174, 3, 62, 170,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2192, 3, 65, 176,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2210, 3, 68, 182,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2228, 3, 71, 188,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2246, 3, 74, 194,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2264, 3, 86, 212,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2282, 3, 89, 218,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2300, 3, 92, 224,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2318, 3, 95, 230,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2336, 3, 98, 236,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2354, 3, 101, 242,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2372, 3, 104, 248,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2390, 3, 107, 254,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2408, 3, 110, 260,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2426, 3, 113, 266,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2444, 3, 116, 272,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2462, 3, 134, 298,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2492, 3, 140, 308,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2522, 3, 146, 318,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2552, 3, 152, 328,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2582, 3, 158, 338,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2612, 3, 164, 348,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2642, 3, 170, 358,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2672, 3, 176, 368,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2702, 3, 182, 378,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2732, 3, 188, 388,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2762, 3, 212, 418,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2792, 3, 218, 428,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2822, 3, 224, 438,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2852, 3, 230, 448,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2882, 3, 236, 458,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2912, 3, 242, 468,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2942, 3, 248, 478,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2972, 3, 254, 488,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3002, 3, 260, 498,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3032, 3, 266, 508,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3062, 3, 298, 548,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3107, 3, 308, 563,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3152, 3, 318, 578,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3197, 3, 328, 593,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3242, 3, 338, 608,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3287, 3, 348, 623,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3332, 3, 358, 638,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3377, 3, 368, 653,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3422, 3, 378, 668,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3467, 3, 418, 713,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3512, 3, 428, 728,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3557, 3, 438, 743,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3602, 3, 448, 758,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3647, 3, 458, 773,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3692, 3, 468, 788,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3737, 3, 478, 803,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3782, 3, 488, 818,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3827, 3, 498, 833,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3872, 3, 548, 890,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3935, 3, 563, 911,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3998, 3, 578, 932,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4061, 3, 593, 953,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4124, 3, 608, 974,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4187, 3, 623, 995,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4250, 3, 638,
                                                                       1016, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4313, 3, 653,
                                                                       1037, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4376, 3, 713,
                                                                       1100, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4439, 3, 728,
                                                                       1121, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4502, 3, 743,
                                                                       1142, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4565, 3, 758,
                                                                       1163, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4628, 3, 773,
                                                                       1184, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4691, 3, 788,
                                                                       1205, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4754, 3, 803,
                                                                       1226, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4817, 3, 818,
                                                                       1247, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4880, 3, 890,
                                                                       1324, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4964, 3, 911,
                                                                       1352, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5048, 3, 932,
                                                                       1380, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5132, 3, 953,
                                                                       1408, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5216, 3, 974,
                                                                       1436, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5300, 3, 995,
                                                                       1464, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5384, 3, 1016,
                                                                       1492, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5468, 3, 1100,
                                                                       1576, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5552, 3, 1121,
                                                                       1604, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5636, 3, 1142,
                                                                       1632, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5720, 3, 1163,
                                                                       1660, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5804, 3, 1184,
                                                                       1688, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5888, 3, 1205,
                                                                       1716, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5972, 3, 1226,
                                                                       1744, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6056, 3, 7, 8,
                                                                       1772, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6062, 3, 8, 9,
                                                                       1775, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6068, 3, 9, 10,
                                                                       1778, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6074, 3, 10, 11,
                                                                       1781, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6080, 3, 11, 12,
                                                                       1784, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6086, 3, 12, 13,
                                                                       1787, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6092, 3, 13, 14,
                                                                       1790, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6098, 3, 14, 15,
                                                                       1793, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6104, 3, 15, 16,
                                                                       1796, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6110, 3, 16, 17,
                                                                       1799, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6116, 3, 17, 18,
                                                                       1802, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6122, 3, 18, 19,
                                                                       1805, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6128, 3, 19, 20,
                                                                       1808, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6134, 3, 23, 24,
                                                                       1811, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6140, 3, 24, 25,
                                                                       1814, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6146, 3, 25, 26,
                                                                       1817, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6152, 3, 26, 27,
                                                                       1820, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6158, 3, 27, 28,
                                                                       1823, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6164, 3, 28, 29,
                                                                       1826, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6170, 3, 29, 30,
                                                                       1829, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6176, 3, 30, 31,
                                                                       1832, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6182, 3, 31, 32,
                                                                       1835, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6188, 3, 32, 33,
                                                                       1838, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6194, 3, 33, 34,
                                                                       1841, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6200, 3, 34, 35,
                                                                       1844, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6206, 3, 35, 36,
                                                                       1847, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6212, 0, 3, 6056,
                                                                       1772, 6062, 1850, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6230, 0, 3, 6062,
                                                                       1775, 6068, 1859, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6248, 0, 3, 6068,
                                                                       1778, 6074, 1868, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6266, 0, 3, 6074,
                                                                       1781, 6080, 1877, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6284, 0, 3, 6080,
                                                                       1784, 6086, 1886, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6302, 0, 3, 6086,
                                                                       1787, 6092, 1895, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6320, 0, 3, 6092,
                                                                       1790, 6098, 1904, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6338, 0, 3, 6098,
                                                                       1793, 6104, 1913, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6356, 0, 3, 6104,
                                                                       1796, 6110, 1922, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6374, 0, 3, 6110,
                                                                       1799, 6116, 1931, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6392, 0, 3, 6116,
                                                                       1802, 6122, 1940, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6410, 0, 3, 6122,
                                                                       1805, 6128, 1949, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6428, 0, 3, 6134,
                                                                       1811, 6140, 1958, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6446, 0, 3, 6140,
                                                                       1814, 6146, 1967, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6464, 0, 3, 6146,
                                                                       1817, 6152, 1976, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6482, 0, 3, 6152,
                                                                       1820, 6158, 1985, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6500, 0, 3, 6158,
                                                                       1823, 6164, 1994, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6518, 0, 3, 6164,
                                                                       1826, 6170, 2003, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6536, 0, 3, 6170,
                                                                       1829, 6176, 2012, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6554, 0, 3, 6176,
                                                                       1832, 6182, 2021, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6572, 0, 3, 6182,
                                                                       1835, 6188, 2030, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6590, 0, 3, 6188,
                                                                       1838, 6194, 2039, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6608, 0, 3, 6194,
                                                                       1841, 6200, 2048, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6626, 0, 3, 6200,
                                                                       1844, 6206, 2057, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6644, 0, 3, 6212,
                                                                       1850, 6230, 122, 128,
                                                                       2066, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6680, 0, 3, 6230,
                                                                       1859, 6248, 128, 134,
                                                                       2084, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6716, 0, 3, 6248,
                                                                       1868, 6266, 134, 140,
                                                                       2102, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6752, 0, 3, 6266,
                                                                       1877, 6284, 140, 146,
                                                                       2120, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6788, 0, 3, 6284,
                                                                       1886, 6302, 146, 152,
                                                                       2138, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6824, 0, 3, 6302,
                                                                       1895, 6320, 152, 158,
                                                                       2156, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6860, 0, 3, 6320,
                                                                       1904, 6338, 158, 164,
                                                                       2174, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6896, 0, 3, 6338,
                                                                       1913, 6356, 164, 170,
                                                                       2192, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6932, 0, 3, 6356,
                                                                       1922, 6374, 170, 176,
                                                                       2210, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6968, 0, 3, 6374,
                                                                       1931, 6392, 176, 182,
                                                                       2228, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7004, 0, 3, 6392,
                                                                       1940, 6410, 182, 188,
                                                                       2246, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7040, 0, 3, 6428,
                                                                       1958, 6446, 200, 206,
                                                                       2264, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7076, 0, 3, 6446,
                                                                       1967, 6464, 206, 212,
                                                                       2282, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7112, 0, 3, 6464,
                                                                       1976, 6482, 212, 218,
                                                                       2300, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7148, 0, 3, 6482,
                                                                       1985, 6500, 218, 224,
                                                                       2318, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7184, 0, 3, 6500,
                                                                       1994, 6518, 224, 230,
                                                                       2336, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7220, 0, 3, 6518,
                                                                       2003, 6536, 230, 236,
                                                                       2354, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7256, 0, 3, 6536,
                                                                       2012, 6554, 236, 242,
                                                                       2372, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7292, 0, 3, 6554,
                                                                       2021, 6572, 242, 248,
                                                                       2390, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7328, 0, 3, 6572,
                                                                       2030, 6590, 248, 254,
                                                                       2408, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7364, 0, 3, 6590,
                                                                       2039, 6608, 254, 260,
                                                                       2426, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7400, 0, 3, 6608,
                                                                       2048, 6626, 260, 266,
                                                                       2444, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7436, 0, 3, 6644,
                                                                       2066, 6680, 278, 288,
                                                                       2462, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7496, 0, 3, 6680,
                                                                       2084, 6716, 288, 298,
                                                                       2492, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7556, 0, 3, 6716,
                                                                       2102, 6752, 298, 308,
                                                                       2522, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7616, 0, 3, 6752,
                                                                       2120, 6788, 308, 318,
                                                                       2552, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7676, 0, 3, 6788,
                                                                       2138, 6824, 318, 328,
                                                                       2582, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7736, 0, 3, 6824,
                                                                       2156, 6860, 328, 338,
                                                                       2612, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7796, 0, 3, 6860,
                                                                       2174, 6896, 338, 348,
                                                                       2642, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7856, 0, 3, 6896,
                                                                       2192, 6932, 348, 358,
                                                                       2672, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7916, 0, 3, 6932,
                                                                       2210, 6968, 358, 368,
                                                                       2702, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7976, 0, 3, 6968,
                                                                       2228, 7004, 368, 378,
                                                                       2732, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8036, 0, 3, 7040,
                                                                       2264, 7076, 398, 408,
                                                                       2762, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8096, 0, 3, 7076,
                                                                       2282, 7112, 408, 418,
                                                                       2792, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8156, 0, 3, 7112,
                                                                       2300, 7148, 418, 428,
                                                                       2822, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8216, 0, 3, 7148,
                                                                       2318, 7184, 428, 438,
                                                                       2852, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8276, 0, 3, 7184,
                                                                       2336, 7220, 438, 448,
                                                                       2882, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8336, 0, 3, 7220,
                                                                       2354, 7256, 448, 458,
                                                                       2912, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8396, 0, 3, 7256,
                                                                       2372, 7292, 458, 468,
                                                                       2942, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8456, 0, 3, 7292,
                                                                       2390, 7328, 468, 478,
                                                                       2972, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8516, 0, 3, 7328,
                                                                       2408, 7364, 478, 488,
                                                                       3002, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8576, 0, 3, 7364,
                                                                       2426, 7400, 488, 498,
                                                                       3032, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8636, 0, 3, 7436,
                                                                       2462, 7496, 518, 533,
                                                                       3062, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8726, 0, 3, 7496,
                                                                       2492, 7556, 533, 548,
                                                                       3107, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8816, 0, 3, 7556,
                                                                       2522, 7616, 548, 563,
                                                                       3152, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8906, 0, 3, 7616,
                                                                       2552, 7676, 563, 578,
                                                                       3197, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8996, 0, 3, 7676,
                                                                       2582, 7736, 578, 593,
                                                                       3242, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9086, 0, 3, 7736,
                                                                       2612, 7796, 593, 608,
                                                                       3287, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9176, 0, 3, 7796,
                                                                       2642, 7856, 608, 623,
                                                                       3332, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9266, 0, 3, 7856,
                                                                       2672, 7916, 623, 638,
                                                                       3377, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9356, 0, 3, 7916,
                                                                       2702, 7976, 638, 653,
                                                                       3422, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9446, 0, 3, 8036,
                                                                       2762, 8096, 683, 698,
                                                                       3467, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9536, 0, 3, 8096,
                                                                       2792, 8156, 698, 713,
                                                                       3512, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9626, 0, 3, 8156,
                                                                       2822, 8216, 713, 728,
                                                                       3557, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9716, 0, 3, 8216,
                                                                       2852, 8276, 728, 743,
                                                                       3602, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9806, 0, 3, 8276,
                                                                       2882, 8336, 743, 758,
                                                                       3647, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9896, 0, 3, 8336,
                                                                       2912, 8396, 758, 773,
                                                                       3692, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9986, 0, 3, 8396,
                                                                       2942, 8456, 773, 788,
                                                                       3737, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10076, 0, 3, 8456,
                                                                       2972, 8516, 788, 803,
                                                                       3782, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10166, 0, 3, 8516,
                                                                       3002, 8576, 803, 818,
                                                                       3827, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10256, 0, 3, 8636,
                                                                       3062, 8726, 848, 869,
                                                                       3872, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10382, 0, 3, 8726,
                                                                       3107, 8816, 869, 890,
                                                                       3935, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10508, 0, 3, 8816,
                                                                       3152, 8906, 890, 911,
                                                                       3998, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10634, 0, 3, 8906,
                                                                       3197, 8996, 911, 932,
                                                                       4061, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10760, 0, 3, 8996,
                                                                       3242, 9086, 932, 953,
                                                                       4124, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10886, 0, 3, 9086,
                                                                       3287, 9176, 953, 974,
                                                                       4187, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11012, 0, 3, 9176,
                                                                       3332, 9266, 974, 995,
                                                                       4250, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11138, 0, 3, 9266,
                                                                       3377, 9356, 995, 1016,
                                                                       4313, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11264, 0, 3, 9446,
                                                                       3467, 9536, 1058, 1079,
                                                                       4376, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11390, 0, 3, 9536,
                                                                       3512, 9626, 1079, 1100,
                                                                       4439, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11516, 0, 3, 9626,
                                                                       3557, 9716, 1100, 1121,
                                                                       4502, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11642, 0, 3, 9716,
                                                                       3602, 9806, 1121, 1142,
                                                                       4565, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11768, 0, 3, 9806,
                                                                       3647, 9896, 1142, 1163,
                                                                       4628, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11894, 0, 3, 9896,
                                                                       3692, 9986, 1163, 1184,
                                                                       4691, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 12020, 0, 3, 9986,
                                                                       3737, 10076, 1184, 1205,
                                                                       4754, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 12146, 0, 3,
                                                                       10076, 3782, 10166, 1205,
                                                                       1226, 4817, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 12272, 0, 3,
                                                                       10256, 3872, 10382, 1268,
                                                                       1296, 4880, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 12440, 0, 3,
                                                                       10382, 3935, 10508, 1296,
                                                                       1324, 4964, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 12608, 0, 3,
                                                                       10508, 3998, 10634, 1324,
                                                                       1352, 5048, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 12776, 0, 3,
                                                                       10634, 4061, 10760, 1352,
                                                                       1380, 5132, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 12944, 0, 3,
                                                                       10760, 4124, 10886, 1380,
                                                                       1408, 5216, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 13112, 0, 3,
                                                                       10886, 4187, 11012, 1408,
                                                                       1436, 5300, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 13280, 0, 3,
                                                                       11012, 4250, 11138, 1436,
                                                                       1464, 5384, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 13448, 0, 3,
                                                                       11264, 4376, 11390, 1520,
                                                                       1548, 5468, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 13616, 0, 3,
                                                                       11390, 4439, 11516, 1548,
                                                                       1576, 5552, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 13784, 0, 3,
                                                                       11516, 4502, 11642, 1576,
                                                                       1604, 5636, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 13952, 0, 3,
                                                                       11642, 4565, 11768, 1604,
                                                                       1632, 5720, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 14120, 0, 3,
                                                                       11768, 4628, 11894, 1632,
                                                                       1660, 5804, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 14288, 0, 3,
                                                                       11894, 4691, 12020, 1660,
                                                                       1688, 5888, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 14456, 0, 3,
                                                                       12020, 4754, 12146, 1688,
                                                                       1716, 5972, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14624, 3, 1772,
                                                                       1775, 6068, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14634, 3, 1775,
                                                                       1778, 6074, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14644, 3, 1778,
                                                                       1781, 6080, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14654, 3, 1781,
                                                                       1784, 6086, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14664, 3, 1784,
                                                                       1787, 6092, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14674, 3, 1787,
                                                                       1790, 6098, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14684, 3, 1790,
                                                                       1793, 6104, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14694, 3, 1793,
                                                                       1796, 6110, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14704, 3, 1796,
                                                                       1799, 6116, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14714, 3, 1799,
                                                                       1802, 6122, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14724, 3, 1802,
                                                                       1805, 6128, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14734, 3, 1811,
                                                                       1814, 6146, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14744, 3, 1814,
                                                                       1817, 6152, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14754, 3, 1817,
                                                                       1820, 6158, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14764, 3, 1820,
                                                                       1823, 6164, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14774, 3, 1823,
                                                                       1826, 6170, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14784, 3, 1826,
                                                                       1829, 6176, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14794, 3, 1829,
                                                                       1832, 6182, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14804, 3, 1832,
                                                                       1835, 6188, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14814, 3, 1835,
                                                                       1838, 6194, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14824, 3, 1838,
                                                                       1841, 6200, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14834, 3, 1841,
                                                                       1844, 6206, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 14844, 0, 3,
                                                                       14624, 6068, 14634, 6248,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 14874, 0, 3,
                                                                       14634, 6074, 14644, 6266,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 14904, 0, 3,
                                                                       14644, 6080, 14654, 6284,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 14934, 0, 3,
                                                                       14654, 6086, 14664, 6302,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 14964, 0, 3,
                                                                       14664, 6092, 14674, 6320,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 14994, 0, 3,
                                                                       14674, 6098, 14684, 6338,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 15024, 0, 3,
                                                                       14684, 6104, 14694, 6356,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 15054, 0, 3,
                                                                       14694, 6110, 14704, 6374,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 15084, 0, 3,
                                                                       14704, 6116, 14714, 6392,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 15114, 0, 3,
                                                                       14714, 6122, 14724, 6410,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 15144, 0, 3,
                                                                       14734, 6146, 14744, 6464,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 15174, 0, 3,
                                                                       14744, 6152, 14754, 6482,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 15204, 0, 3,
                                                                       14754, 6158, 14764, 6500,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 15234, 0, 3,
                                                                       14764, 6164, 14774, 6518,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 15264, 0, 3,
                                                                       14774, 6170, 14784, 6536,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 15294, 0, 3,
                                                                       14784, 6176, 14794, 6554,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 15324, 0, 3,
                                                                       14794, 6182, 14804, 6572,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 15354, 0, 3,
                                                                       14804, 6188, 14814, 6590,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 15384, 0, 3,
                                                                       14814, 6194, 14824, 6608,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 15414, 0, 3,
                                                                       14824, 6200, 14834, 6626,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 15444, 0, 3,
                                                                       14844, 6248, 14874, 2066,
                                                                       2084, 6716, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 15504, 0, 3,
                                                                       14874, 6266, 14904, 2084,
                                                                       2102, 6752, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 15564, 0, 3,
                                                                       14904, 6284, 14934, 2102,
                                                                       2120, 6788, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 15624, 0, 3,
                                                                       14934, 6302, 14964, 2120,
                                                                       2138, 6824, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 15684, 0, 3,
                                                                       14964, 6320, 14994, 2138,
                                                                       2156, 6860, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 15744, 0, 3,
                                                                       14994, 6338, 15024, 2156,
                                                                       2174, 6896, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 15804, 0, 3,
                                                                       15024, 6356, 15054, 2174,
                                                                       2192, 6932, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 15864, 0, 3,
                                                                       15054, 6374, 15084, 2192,
                                                                       2210, 6968, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 15924, 0, 3,
                                                                       15084, 6392, 15114, 2210,
                                                                       2228, 7004, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 15984, 0, 3,
                                                                       15144, 6464, 15174, 2264,
                                                                       2282, 7112, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 16044, 0, 3,
                                                                       15174, 6482, 15204, 2282,
                                                                       2300, 7148, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 16104, 0, 3,
                                                                       15204, 6500, 15234, 2300,
                                                                       2318, 7184, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 16164, 0, 3,
                                                                       15234, 6518, 15264, 2318,
                                                                       2336, 7220, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 16224, 0, 3,
                                                                       15264, 6536, 15294, 2336,
                                                                       2354, 7256, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 16284, 0, 3,
                                                                       15294, 6554, 15324, 2354,
                                                                       2372, 7292, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 16344, 0, 3,
                                                                       15324, 6572, 15354, 2372,
                                                                       2390, 7328, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 16404, 0, 3,
                                                                       15354, 6590, 15384, 2390,
                                                                       2408, 7364, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 16464, 0, 3,
                                                                       15384, 6608, 15414, 2408,
                                                                       2426, 7400, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 16524, 0, 3,
                                                                       15444, 6716, 15504, 2462,
                                                                       2492, 7556, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 16624, 0, 3,
                                                                       15504, 6752, 15564, 2492,
                                                                       2522, 7616, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 16724, 0, 3,
                                                                       15564, 6788, 15624, 2522,
                                                                       2552, 7676, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 16824, 0, 3,
                                                                       15624, 6824, 15684, 2552,
                                                                       2582, 7736, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 16924, 0, 3,
                                                                       15684, 6860, 15744, 2582,
                                                                       2612, 7796, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 17024, 0, 3,
                                                                       15744, 6896, 15804, 2612,
                                                                       2642, 7856, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 17124, 0, 3,
                                                                       15804, 6932, 15864, 2642,
                                                                       2672, 7916, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 17224, 0, 3,
                                                                       15864, 6968, 15924, 2672,
                                                                       2702, 7976, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 17324, 0, 3,
                                                                       15984, 7112, 16044, 2762,
                                                                       2792, 8156, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 17424, 0, 3,
                                                                       16044, 7148, 16104, 2792,
                                                                       2822, 8216, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 17524, 0, 3,
                                                                       16104, 7184, 16164, 2822,
                                                                       2852, 8276, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 17624, 0, 3,
                                                                       16164, 7220, 16224, 2852,
                                                                       2882, 8336, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 17724, 0, 3,
                                                                       16224, 7256, 16284, 2882,
                                                                       2912, 8396, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 17824, 0, 3,
                                                                       16284, 7292, 16344, 2912,
                                                                       2942, 8456, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 17924, 0, 3,
                                                                       16344, 7328, 16404, 2942,
                                                                       2972, 8516, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 18024, 0, 3,
                                                                       16404, 7364, 16464, 2972,
                                                                       3002, 8576, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 18124, 0, 3,
                                                                       16524, 7556, 16624, 3062,
                                                                       3107, 8816, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 18274, 0, 3,
                                                                       16624, 7616, 16724, 3107,
                                                                       3152, 8906, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 18424, 0, 3,
                                                                       16724, 7676, 16824, 3152,
                                                                       3197, 8996, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 18574, 0, 3,
                                                                       16824, 7736, 16924, 3197,
                                                                       3242, 9086, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 18724, 0, 3,
                                                                       16924, 7796, 17024, 3242,
                                                                       3287, 9176, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 18874, 0, 3,
                                                                       17024, 7856, 17124, 3287,
                                                                       3332, 9266, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 19024, 0, 3,
                                                                       17124, 7916, 17224, 3332,
                                                                       3377, 9356, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 19174, 0, 3,
                                                                       17324, 8156, 17424, 3467,
                                                                       3512, 9626, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 19324, 0, 3,
                                                                       17424, 8216, 17524, 3512,
                                                                       3557, 9716, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 19474, 0, 3,
                                                                       17524, 8276, 17624, 3557,
                                                                       3602, 9806, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 19624, 0, 3,
                                                                       17624, 8336, 17724, 3602,
                                                                       3647, 9896, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 19774, 0, 3,
                                                                       17724, 8396, 17824, 3647,
                                                                       3692, 9986, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 19924, 0, 3,
                                                                       17824, 8456, 17924, 3692,
                                                                       3737, 10076, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 20074, 0, 3,
                                                                       17924, 8516, 18024, 3737,
                                                                       3782, 10166, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 20224, 0, 3,
                                                                       18124, 8816, 18274, 3872,
                                                                       3935, 10508, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 20434, 0, 3,
                                                                       18274, 8906, 18424, 3935,
                                                                       3998, 10634, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 20644, 0, 3,
                                                                       18424, 8996, 18574, 3998,
                                                                       4061, 10760, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 20854, 0, 3,
                                                                       18574, 9086, 18724, 4061,
                                                                       4124, 10886, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 21064, 0, 3,
                                                                       18724, 9176, 18874, 4124,
                                                                       4187, 11012, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 21274, 0, 3,
                                                                       18874, 9266, 19024, 4187,
                                                                       4250, 11138, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 21484, 0, 3,
                                                                       19174, 9626, 19324, 4376,
                                                                       4439, 11516, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 21694, 0, 3,
                                                                       19324, 9716, 19474, 4439,
                                                                       4502, 11642, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 21904, 0, 3,
                                                                       19474, 9806, 19624, 4502,
                                                                       4565, 11768, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 22114, 0, 3,
                                                                       19624, 9896, 19774, 4565,
                                                                       4628, 11894, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 22324, 0, 3,
                                                                       19774, 9986, 19924, 4628,
                                                                       4691, 12020, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 22534, 0, 3,
                                                                       19924, 10076, 20074, 4691,
                                                                       4754, 12146, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 22744, 0, 3,
                                                                       20224, 10508, 20434, 4880,
                                                                       4964, 12608, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 23024, 0, 3,
                                                                       20434, 10634, 20644, 4964,
                                                                       5048, 12776, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 23304, 0, 3,
                                                                       20644, 10760, 20854, 5048,
                                                                       5132, 12944, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 23584, 0, 3,
                                                                       20854, 10886, 21064, 5132,
                                                                       5216, 13112, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 23864, 0, 3,
                                                                       21064, 11012, 21274, 5216,
                                                                       5300, 13280, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 24144, 0, 3,
                                                                       21484, 11516, 21694, 5468,
                                                                       5552, 13784, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 24424, 0, 3,
                                                                       21694, 11642, 21904, 5552,
                                                                       5636, 13952, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 24704, 0, 3,
                                                                       21904, 11768, 22114, 5636,
                                                                       5720, 14120, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 24984, 0, 3,
                                                                       22114, 11894, 22324, 5720,
                                                                       5804, 14288, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 25264, 0, 3,
                                                                       22324, 12020, 22534, 5804,
                                                                       5888, 14456, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 25544, 3, 6056,
                                                                       6062, 14624, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 25559, 3, 6062,
                                                                       6068, 14634, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 25574, 3, 6068,
                                                                       6074, 14644, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 25589, 3, 6074,
                                                                       6080, 14654, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 25604, 3, 6080,
                                                                       6086, 14664, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 25619, 3, 6086,
                                                                       6092, 14674, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 25634, 3, 6092,
                                                                       6098, 14684, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 25649, 3, 6098,
                                                                       6104, 14694, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 25664, 3, 6104,
                                                                       6110, 14704, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 25679, 3, 6110,
                                                                       6116, 14714, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 25694, 3, 6116,
                                                                       6122, 14724, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 25709, 3, 6134,
                                                                       6140, 14734, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 25724, 3, 6140,
                                                                       6146, 14744, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 25739, 3, 6146,
                                                                       6152, 14754, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 25754, 3, 6152,
                                                                       6158, 14764, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 25769, 3, 6158,
                                                                       6164, 14774, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 25784, 3, 6164,
                                                                       6170, 14784, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 25799, 3, 6170,
                                                                       6176, 14794, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 25814, 3, 6176,
                                                                       6182, 14804, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 25829, 3, 6182,
                                                                       6188, 14814, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 25844, 3, 6188,
                                                                       6194, 14824, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 25859, 3, 6194,
                                                                       6200, 14834, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 25874, 0, 3,
                                                                       25544, 14624, 25559, 6212,
                                                                       6230, 14844, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 25919, 0, 3,
                                                                       25559, 14634, 25574, 6230,
                                                                       6248, 14874, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 25964, 0, 3,
                                                                       25574, 14644, 25589, 6248,
                                                                       6266, 14904, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 26009, 0, 3,
                                                                       25589, 14654, 25604, 6266,
                                                                       6284, 14934, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 26054, 0, 3,
                                                                       25604, 14664, 25619, 6284,
                                                                       6302, 14964, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 26099, 0, 3,
                                                                       25619, 14674, 25634, 6302,
                                                                       6320, 14994, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 26144, 0, 3,
                                                                       25634, 14684, 25649, 6320,
                                                                       6338, 15024, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 26189, 0, 3,
                                                                       25649, 14694, 25664, 6338,
                                                                       6356, 15054, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 26234, 0, 3,
                                                                       25664, 14704, 25679, 6356,
                                                                       6374, 15084, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 26279, 0, 3,
                                                                       25679, 14714, 25694, 6374,
                                                                       6392, 15114, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 26324, 0, 3,
                                                                       25709, 14734, 25724, 6428,
                                                                       6446, 15144, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 26369, 0, 3,
                                                                       25724, 14744, 25739, 6446,
                                                                       6464, 15174, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 26414, 0, 3,
                                                                       25739, 14754, 25754, 6464,
                                                                       6482, 15204, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 26459, 0, 3,
                                                                       25754, 14764, 25769, 6482,
                                                                       6500, 15234, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 26504, 0, 3,
                                                                       25769, 14774, 25784, 6500,
                                                                       6518, 15264, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 26549, 0, 3,
                                                                       25784, 14784, 25799, 6518,
                                                                       6536, 15294, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 26594, 0, 3,
                                                                       25799, 14794, 25814, 6536,
                                                                       6554, 15324, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 26639, 0, 3,
                                                                       25814, 14804, 25829, 6554,
                                                                       6572, 15354, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 26684, 0, 3,
                                                                       25829, 14814, 25844, 6572,
                                                                       6590, 15384, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 26729, 0, 3,
                                                                       25844, 14824, 25859, 6590,
                                                                       6608, 15414, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 26774, 0, 3,
                                                                       25874, 14844, 25919, 6644,
                                                                       6680, 15444, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 26864, 0, 3,
                                                                       25919, 14874, 25964, 6680,
                                                                       6716, 15504, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 26954, 0, 3,
                                                                       25964, 14904, 26009, 6716,
                                                                       6752, 15564, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 27044, 0, 3,
                                                                       26009, 14934, 26054, 6752,
                                                                       6788, 15624, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 27134, 0, 3,
                                                                       26054, 14964, 26099, 6788,
                                                                       6824, 15684, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 27224, 0, 3,
                                                                       26099, 14994, 26144, 6824,
                                                                       6860, 15744, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 27314, 0, 3,
                                                                       26144, 15024, 26189, 6860,
                                                                       6896, 15804, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 27404, 0, 3,
                                                                       26189, 15054, 26234, 6896,
                                                                       6932, 15864, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 27494, 0, 3,
                                                                       26234, 15084, 26279, 6932,
                                                                       6968, 15924, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 27584, 0, 3,
                                                                       26324, 15144, 26369, 7040,
                                                                       7076, 15984, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 27674, 0, 3,
                                                                       26369, 15174, 26414, 7076,
                                                                       7112, 16044, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 27764, 0, 3,
                                                                       26414, 15204, 26459, 7112,
                                                                       7148, 16104, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 27854, 0, 3,
                                                                       26459, 15234, 26504, 7148,
                                                                       7184, 16164, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 27944, 0, 3,
                                                                       26504, 15264, 26549, 7184,
                                                                       7220, 16224, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 28034, 0, 3,
                                                                       26549, 15294, 26594, 7220,
                                                                       7256, 16284, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 28124, 0, 3,
                                                                       26594, 15324, 26639, 7256,
                                                                       7292, 16344, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 28214, 0, 3,
                                                                       26639, 15354, 26684, 7292,
                                                                       7328, 16404, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 28304, 0, 3,
                                                                       26684, 15384, 26729, 7328,
                                                                       7364, 16464, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 28394, 0, 3,
                                                                       26774, 15444, 26864, 7436,
                                                                       7496, 16524, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 28544, 0, 3,
                                                                       26864, 15504, 26954, 7496,
                                                                       7556, 16624, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 28694, 0, 3,
                                                                       26954, 15564, 27044, 7556,
                                                                       7616, 16724, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 28844, 0, 3,
                                                                       27044, 15624, 27134, 7616,
                                                                       7676, 16824, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 28994, 0, 3,
                                                                       27134, 15684, 27224, 7676,
                                                                       7736, 16924, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 29144, 0, 3,
                                                                       27224, 15744, 27314, 7736,
                                                                       7796, 17024, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 29294, 0, 3,
                                                                       27314, 15804, 27404, 7796,
                                                                       7856, 17124, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 29444, 0, 3,
                                                                       27404, 15864, 27494, 7856,
                                                                       7916, 17224, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 29594, 0, 3,
                                                                       27584, 15984, 27674, 8036,
                                                                       8096, 17324, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 29744, 0, 3,
                                                                       27674, 16044, 27764, 8096,
                                                                       8156, 17424, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 29894, 0, 3,
                                                                       27764, 16104, 27854, 8156,
                                                                       8216, 17524, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 30044, 0, 3,
                                                                       27854, 16164, 27944, 8216,
                                                                       8276, 17624, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 30194, 0, 3,
                                                                       27944, 16224, 28034, 8276,
                                                                       8336, 17724, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 30344, 0, 3,
                                                                       28034, 16284, 28124, 8336,
                                                                       8396, 17824, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 30494, 0, 3,
                                                                       28124, 16344, 28214, 8396,
                                                                       8456, 17924, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 30644, 0, 3,
                                                                       28214, 16404, 28304, 8456,
                                                                       8516, 18024, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 30794, 0, 3,
                                                                       28394, 16524, 28544, 8636,
                                                                       8726, 18124, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 31019, 0, 3,
                                                                       28544, 16624, 28694, 8726,
                                                                       8816, 18274, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 31244, 0, 3,
                                                                       28694, 16724, 28844, 8816,
                                                                       8906, 18424, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 31469, 0, 3,
                                                                       28844, 16824, 28994, 8906,
                                                                       8996, 18574, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 31694, 0, 3,
                                                                       28994, 16924, 29144, 8996,
                                                                       9086, 18724, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 31919, 0, 3,
                                                                       29144, 17024, 29294, 9086,
                                                                       9176, 18874, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 32144, 0, 3,
                                                                       29294, 17124, 29444, 9176,
                                                                       9266, 19024, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 32369, 0, 3,
                                                                       29594, 17324, 29744, 9446,
                                                                       9536, 19174, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 32594, 0, 3,
                                                                       29744, 17424, 29894, 9536,
                                                                       9626, 19324, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 32819, 0, 3,
                                                                       29894, 17524, 30044, 9626,
                                                                       9716, 19474, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 33044, 0, 3,
                                                                       30044, 17624, 30194, 9716,
                                                                       9806, 19624, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 33269, 0, 3,
                                                                       30194, 17724, 30344, 9806,
                                                                       9896, 19774, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 33494, 0, 3,
                                                                       30344, 17824, 30494, 9896,
                                                                       9986, 19924, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 33719, 0, 3,
                                                                       30494, 17924, 30644, 9986,
                                                                       10076, 20074, ncols,
                                                                       gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 33944, 0, 3,
                                                                       30794, 18124, 31019,
                                                                       10256, 10382, 20224,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 34259, 0, 3,
                                                                       31019, 18274, 31244,
                                                                       10382, 10508, 20434,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 34574, 0, 3,
                                                                       31244, 18424, 31469,
                                                                       10508, 10634, 20644,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 34889, 0, 3,
                                                                       31469, 18574, 31694,
                                                                       10634, 10760, 20854,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 35204, 0, 3,
                                                                       31694, 18724, 31919,
                                                                       10760, 10886, 21064,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 35519, 0, 3,
                                                                       31919, 18874, 32144,
                                                                       10886, 11012, 21274,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 35834, 0, 3,
                                                                       32369, 19174, 32594,
                                                                       11264, 11390, 21484,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 36149, 0, 3,
                                                                       32594, 19324, 32819,
                                                                       11390, 11516, 21694,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 36464, 0, 3,
                                                                       32819, 19474, 33044,
                                                                       11516, 11642, 21904,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 36779, 0, 3,
                                                                       33044, 19624, 33269,
                                                                       11642, 11768, 22114,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 37094, 0, 3,
                                                                       33269, 19774, 33494,
                                                                       11768, 11894, 22324,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 37409, 0, 3,
                                                                       33494, 19924, 33719,
                                                                       11894, 12020, 22534,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 37724, 0, 3,
                                                                       33944, 20224, 34259,
                                                                       12272, 12440, 22744,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 38144, 0, 3,
                                                                       34259, 20434, 34574,
                                                                       12440, 12608, 23024,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 38564, 0, 3,
                                                                       34574, 20644, 34889,
                                                                       12608, 12776, 23304,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 38984, 0, 3,
                                                                       34889, 20854, 35204,
                                                                       12776, 12944, 23584,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 39404, 0, 3,
                                                                       35204, 21064, 35519,
                                                                       12944, 13112, 23864,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 39824, 0, 3,
                                                                       35834, 21484, 36149,
                                                                       13448, 13616, 24144,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 40244, 0, 3,
                                                                       36149, 21694, 36464,
                                                                       13616, 13784, 24424,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 40664, 0, 3,
                                                                       36464, 21904, 36779,
                                                                       13784, 13952, 24704,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 41084, 0, 3,
                                                                       36779, 22114, 37094,
                                                                       13952, 14120, 24984,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 41504, 0, 3,
                                                                       37094, 22324, 37409,
                                                                       14120, 14288, 25264,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 41924, 3, 14624,
                                                                       14634, 25574, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 41945, 3, 14634,
                                                                       14644, 25589, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 41966, 3, 14644,
                                                                       14654, 25604, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 41987, 3, 14654,
                                                                       14664, 25619, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 42008, 3, 14664,
                                                                       14674, 25634, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 42029, 3, 14674,
                                                                       14684, 25649, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 42050, 3, 14684,
                                                                       14694, 25664, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 42071, 3, 14694,
                                                                       14704, 25679, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 42092, 3, 14704,
                                                                       14714, 25694, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 42113, 3, 14734,
                                                                       14744, 25739, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 42134, 3, 14744,
                                                                       14754, 25754, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 42155, 3, 14754,
                                                                       14764, 25769, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 42176, 3, 14764,
                                                                       14774, 25784, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 42197, 3, 14774,
                                                                       14784, 25799, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 42218, 3, 14784,
                                                                       14794, 25814, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 42239, 3, 14794,
                                                                       14804, 25829, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 42260, 3, 14804,
                                                                       14814, 25844, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 42281, 3, 14814,
                                                                       14824, 25859, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 42302, 0, 3,
                                                                       41924, 25574, 41945,
                                                                       14844, 14874, 25964,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 42365, 0, 3,
                                                                       41945, 25589, 41966,
                                                                       14874, 14904, 26009,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 42428, 0, 3,
                                                                       41966, 25604, 41987,
                                                                       14904, 14934, 26054,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 42491, 0, 3,
                                                                       41987, 25619, 42008,
                                                                       14934, 14964, 26099,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 42554, 0, 3,
                                                                       42008, 25634, 42029,
                                                                       14964, 14994, 26144,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 42617, 0, 3,
                                                                       42029, 25649, 42050,
                                                                       14994, 15024, 26189,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 42680, 0, 3,
                                                                       42050, 25664, 42071,
                                                                       15024, 15054, 26234,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 42743, 0, 3,
                                                                       42071, 25679, 42092,
                                                                       15054, 15084, 26279,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 42806, 0, 3,
                                                                       42113, 25739, 42134,
                                                                       15144, 15174, 26414,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 42869, 0, 3,
                                                                       42134, 25754, 42155,
                                                                       15174, 15204, 26459,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 42932, 0, 3,
                                                                       42155, 25769, 42176,
                                                                       15204, 15234, 26504,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 42995, 0, 3,
                                                                       42176, 25784, 42197,
                                                                       15234, 15264, 26549,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 43058, 0, 3,
                                                                       42197, 25799, 42218,
                                                                       15264, 15294, 26594,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 43121, 0, 3,
                                                                       42218, 25814, 42239,
                                                                       15294, 15324, 26639,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 43184, 0, 3,
                                                                       42239, 25829, 42260,
                                                                       15324, 15354, 26684,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 43247, 0, 3,
                                                                       42260, 25844, 42281,
                                                                       15354, 15384, 26729,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 43310, 0, 3,
                                                                       42302, 25964, 42365,
                                                                       15444, 15504, 26954,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 43436, 0, 3,
                                                                       42365, 26009, 42428,
                                                                       15504, 15564, 27044,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 43562, 0, 3,
                                                                       42428, 26054, 42491,
                                                                       15564, 15624, 27134,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 43688, 0, 3,
                                                                       42491, 26099, 42554,
                                                                       15624, 15684, 27224,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 43814, 0, 3,
                                                                       42554, 26144, 42617,
                                                                       15684, 15744, 27314,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 43940, 0, 3,
                                                                       42617, 26189, 42680,
                                                                       15744, 15804, 27404,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 44066, 0, 3,
                                                                       42680, 26234, 42743,
                                                                       15804, 15864, 27494,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 44192, 0, 3,
                                                                       42806, 26414, 42869,
                                                                       15984, 16044, 27764,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 44318, 0, 3,
                                                                       42869, 26459, 42932,
                                                                       16044, 16104, 27854,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 44444, 0, 3,
                                                                       42932, 26504, 42995,
                                                                       16104, 16164, 27944,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 44570, 0, 3,
                                                                       42995, 26549, 43058,
                                                                       16164, 16224, 28034,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 44696, 0, 3,
                                                                       43058, 26594, 43121,
                                                                       16224, 16284, 28124,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 44822, 0, 3,
                                                                       43121, 26639, 43184,
                                                                       16284, 16344, 28214,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 44948, 0, 3,
                                                                       43184, 26684, 43247,
                                                                       16344, 16404, 28304,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 45074, 0, 3,
                                                                       43310, 26954, 43436,
                                                                       16524, 16624, 28694,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 45284, 0, 3,
                                                                       43436, 27044, 43562,
                                                                       16624, 16724, 28844,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 45494, 0, 3,
                                                                       43562, 27134, 43688,
                                                                       16724, 16824, 28994,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 45704, 0, 3,
                                                                       43688, 27224, 43814,
                                                                       16824, 16924, 29144,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 45914, 0, 3,
                                                                       43814, 27314, 43940,
                                                                       16924, 17024, 29294,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 46124, 0, 3,
                                                                       43940, 27404, 44066,
                                                                       17024, 17124, 29444,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 46334, 0, 3,
                                                                       44192, 27764, 44318,
                                                                       17324, 17424, 29894,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 46544, 0, 3,
                                                                       44318, 27854, 44444,
                                                                       17424, 17524, 30044,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 46754, 0, 3,
                                                                       44444, 27944, 44570,
                                                                       17524, 17624, 30194,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 46964, 0, 3,
                                                                       44570, 28034, 44696,
                                                                       17624, 17724, 30344,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 47174, 0, 3,
                                                                       44696, 28124, 44822,
                                                                       17724, 17824, 30494,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 47384, 0, 3,
                                                                       44822, 28214, 44948,
                                                                       17824, 17924, 30644,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 47594, 0, 3,
                                                                       45074, 28694, 45284,
                                                                       18124, 18274, 31244,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 47909, 0, 3,
                                                                       45284, 28844, 45494,
                                                                       18274, 18424, 31469,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 48224, 0, 3,
                                                                       45494, 28994, 45704,
                                                                       18424, 18574, 31694,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 48539, 0, 3,
                                                                       45704, 29144, 45914,
                                                                       18574, 18724, 31919,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 48854, 0, 3,
                                                                       45914, 29294, 46124,
                                                                       18724, 18874, 32144,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 49169, 0, 3,
                                                                       46334, 29894, 46544,
                                                                       19174, 19324, 32819,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 49484, 0, 3,
                                                                       46544, 30044, 46754,
                                                                       19324, 19474, 33044,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 49799, 0, 3,
                                                                       46754, 30194, 46964,
                                                                       19474, 19624, 33269,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 50114, 0, 3,
                                                                       46964, 30344, 47174,
                                                                       19624, 19774, 33494,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 50429, 0, 3,
                                                                       47174, 30494, 47384,
                                                                       19774, 19924, 33719,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 50744, 0, 3,
                                                                       47594, 31244, 47909,
                                                                       20224, 20434, 34574,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 51185, 0, 3,
                                                                       47909, 31469, 48224,
                                                                       20434, 20644, 34889,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 51626, 0, 3,
                                                                       48224, 31694, 48539,
                                                                       20644, 20854, 35204,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 52067, 0, 3,
                                                                       48539, 31919, 48854,
                                                                       20854, 21064, 35519,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 52508, 0, 3,
                                                                       49169, 32819, 49484,
                                                                       21484, 21694, 36464,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 52949, 0, 3,
                                                                       49484, 33044, 49799,
                                                                       21694, 21904, 36779,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 53390, 0, 3,
                                                                       49799, 33269, 50114,
                                                                       21904, 22114, 37094,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 53831, 0, 3,
                                                                       50114, 33494, 50429,
                                                                       22114, 22324, 37409,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 54272, 0, 3,
                                                                       50744, 34574, 51185,
                                                                       22744, 23024, 38564,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 54860, 0, 3,
                                                                       51185, 34889, 51626,
                                                                       23024, 23304, 38984,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 55448, 0, 3,
                                                                       51626, 35204, 52067,
                                                                       23304, 23584, 39404,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 56036, 0, 3,
                                                                       52508, 36464, 52949,
                                                                       24144, 24424, 40664,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 56624, 0, 3,
                                                                       52949, 36779, 53390,
                                                                       24424, 24704, 41084,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 57212, 0, 3,
                                                                       53390, 37094, 53831,
                                                                       24704, 24984, 41504,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 57800, 3, 25544,
                                                                       25559, 41924, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 57828, 3, 25559,
                                                                       25574, 41945, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 57856, 3, 25574,
                                                                       25589, 41966, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 57884, 3, 25589,
                                                                       25604, 41987, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 57912, 3, 25604,
                                                                       25619, 42008, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 57940, 3, 25619,
                                                                       25634, 42029, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 57968, 3, 25634,
                                                                       25649, 42050, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 57996, 3, 25649,
                                                                       25664, 42071, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 58024, 3, 25664,
                                                                       25679, 42092, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 58052, 3, 25709,
                                                                       25724, 42113, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 58080, 3, 25724,
                                                                       25739, 42134, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 58108, 3, 25739,
                                                                       25754, 42155, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 58136, 3, 25754,
                                                                       25769, 42176, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 58164, 3, 25769,
                                                                       25784, 42197, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 58192, 3, 25784,
                                                                       25799, 42218, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 58220, 3, 25799,
                                                                       25814, 42239, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 58248, 3, 25814,
                                                                       25829, 42260, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 58276, 3, 25829,
                                                                       25844, 42281, ncols,
                                                                       gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 58304, 0, 3,
                                                                       57800, 41924, 57828,
                                                                       25874, 25919, 42302,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 58388, 0, 3,
                                                                       57828, 41945, 57856,
                                                                       25919, 25964, 42365,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 58472, 0, 3,
                                                                       57856, 41966, 57884,
                                                                       25964, 26009, 42428,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 58556, 0, 3,
                                                                       57884, 41987, 57912,
                                                                       26009, 26054, 42491,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 58640, 0, 3,
                                                                       57912, 42008, 57940,
                                                                       26054, 26099, 42554,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 58724, 0, 3,
                                                                       57940, 42029, 57968,
                                                                       26099, 26144, 42617,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 58808, 0, 3,
                                                                       57968, 42050, 57996,
                                                                       26144, 26189, 42680,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 58892, 0, 3,
                                                                       57996, 42071, 58024,
                                                                       26189, 26234, 42743,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 58976, 0, 3,
                                                                       58052, 42113, 58080,
                                                                       26324, 26369, 42806,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 59060, 0, 3,
                                                                       58080, 42134, 58108,
                                                                       26369, 26414, 42869,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 59144, 0, 3,
                                                                       58108, 42155, 58136,
                                                                       26414, 26459, 42932,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 59228, 0, 3,
                                                                       58136, 42176, 58164,
                                                                       26459, 26504, 42995,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 59312, 0, 3,
                                                                       58164, 42197, 58192,
                                                                       26504, 26549, 43058,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 59396, 0, 3,
                                                                       58192, 42218, 58220,
                                                                       26549, 26594, 43121,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 59480, 0, 3,
                                                                       58220, 42239, 58248,
                                                                       26594, 26639, 43184,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 59564, 0, 3,
                                                                       58248, 42260, 58276,
                                                                       26639, 26684, 43247,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 59648, 0, 3,
                                                                       58304, 42302, 58388,
                                                                       26774, 26864, 43310,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 59816, 0, 3,
                                                                       58388, 42365, 58472,
                                                                       26864, 26954, 43436,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 59984, 0, 3,
                                                                       58472, 42428, 58556,
                                                                       26954, 27044, 43562,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 60152, 0, 3,
                                                                       58556, 42491, 58640,
                                                                       27044, 27134, 43688,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 60320, 0, 3,
                                                                       58640, 42554, 58724,
                                                                       27134, 27224, 43814,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 60488, 0, 3,
                                                                       58724, 42617, 58808,
                                                                       27224, 27314, 43940,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 60656, 0, 3,
                                                                       58808, 42680, 58892,
                                                                       27314, 27404, 44066,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 60824, 0, 3,
                                                                       58976, 42806, 59060,
                                                                       27584, 27674, 44192,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 60992, 0, 3,
                                                                       59060, 42869, 59144,
                                                                       27674, 27764, 44318,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 61160, 0, 3,
                                                                       59144, 42932, 59228,
                                                                       27764, 27854, 44444,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 61328, 0, 3,
                                                                       59228, 42995, 59312,
                                                                       27854, 27944, 44570,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 61496, 0, 3,
                                                                       59312, 43058, 59396,
                                                                       27944, 28034, 44696,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 61664, 0, 3,
                                                                       59396, 43121, 59480,
                                                                       28034, 28124, 44822,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 61832, 0, 3,
                                                                       59480, 43184, 59564,
                                                                       28124, 28214, 44948,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 62000, 0, 3,
                                                                       59648, 43310, 59816,
                                                                       28394, 28544, 45074,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 62280, 0, 3,
                                                                       59816, 43436, 59984,
                                                                       28544, 28694, 45284,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 62560, 0, 3,
                                                                       59984, 43562, 60152,
                                                                       28694, 28844, 45494,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 62840, 0, 3,
                                                                       60152, 43688, 60320,
                                                                       28844, 28994, 45704,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 63120, 0, 3,
                                                                       60320, 43814, 60488,
                                                                       28994, 29144, 45914,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 63400, 0, 3,
                                                                       60488, 43940, 60656,
                                                                       29144, 29294, 46124,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 63680, 0, 3,
                                                                       60824, 44192, 60992,
                                                                       29594, 29744, 46334,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 63960, 0, 3,
                                                                       60992, 44318, 61160,
                                                                       29744, 29894, 46544,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 64240, 0, 3,
                                                                       61160, 44444, 61328,
                                                                       29894, 30044, 46754,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 64520, 0, 3,
                                                                       61328, 44570, 61496,
                                                                       30044, 30194, 46964,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 64800, 0, 3,
                                                                       61496, 44696, 61664,
                                                                       30194, 30344, 47174,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 65080, 0, 3,
                                                                       61664, 44822, 61832,
                                                                       30344, 30494, 47384,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 65360, 0, 3,
                                                                       62000, 45074, 62280,
                                                                       30794, 31019, 47594,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 65780, 0, 3,
                                                                       62280, 45284, 62560,
                                                                       31019, 31244, 47909,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 66200, 0, 3,
                                                                       62560, 45494, 62840,
                                                                       31244, 31469, 48224,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 66620, 0, 3,
                                                                       62840, 45704, 63120,
                                                                       31469, 31694, 48539,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 67040, 0, 3,
                                                                       63120, 45914, 63400,
                                                                       31694, 31919, 48854,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 67460, 0, 3,
                                                                       63680, 46334, 63960,
                                                                       32369, 32594, 49169,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 67880, 0, 3,
                                                                       63960, 46544, 64240,
                                                                       32594, 32819, 49484,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 68300, 0, 3,
                                                                       64240, 46754, 64520,
                                                                       32819, 33044, 49799,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 68720, 0, 3,
                                                                       64520, 46964, 64800,
                                                                       33044, 33269, 50114,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 69140, 0, 3,
                                                                       64800, 47174, 65080,
                                                                       33269, 33494, 50429,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 69560, 0, 3,
                                                                       65360, 47594, 65780,
                                                                       33944, 34259, 50744,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 70148, 0, 3,
                                                                       65780, 47909, 66200,
                                                                       34259, 34574, 51185,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 70736, 0, 3,
                                                                       66200, 48224, 66620,
                                                                       34574, 34889, 51626,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 71324, 0, 3,
                                                                       66620, 48539, 67040,
                                                                       34889, 35204, 52067,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 71912, 0, 3,
                                                                       67460, 49169, 67880,
                                                                       35834, 36149, 52508,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 72500, 0, 3,
                                                                       67880, 49484, 68300,
                                                                       36149, 36464, 52949,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 73088, 0, 3,
                                                                       68300, 49799, 68720,
                                                                       36464, 36779, 53390,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 73676, 0, 3,
                                                                       68720, 50114, 69140,
                                                                       36779, 37094, 53831,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 74264, 0, 3,
                                                                       69560, 50744, 70148,
                                                                       37724, 38144, 54272,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 75048, 0, 3,
                                                                       70148, 51185, 70736,
                                                                       38144, 38564, 54860,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 75832, 0, 3,
                                                                       70736, 51626, 71324,
                                                                       38564, 38984, 55448,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 76616, 0, 3,
                                                                       71912, 52508, 72500,
                                                                       39824, 40244, 56036,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 77400, 0, 3,
                                                                       72500, 52949, 73088,
                                                                       40244, 40664, 56624,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 78184, 0, 3,
                                                                       73088, 53390, 73676,
                                                                       40664, 41084, 57212,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 78968, 3, 41924,
                                                                       41945, 57856, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 79004, 3, 41945,
                                                                       41966, 57884, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 79040, 3, 41966,
                                                                       41987, 57912, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 79076, 3, 41987,
                                                                       42008, 57940, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 79112, 3, 42008,
                                                                       42029, 57968, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 79148, 3, 42029,
                                                                       42050, 57996, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 79184, 3, 42050,
                                                                       42071, 58024, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 79220, 3, 42113,
                                                                       42134, 58108, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 79256, 3, 42134,
                                                                       42155, 58136, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 79292, 3, 42155,
                                                                       42176, 58164, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 79328, 3, 42176,
                                                                       42197, 58192, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 79364, 3, 42197,
                                                                       42218, 58220, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 79400, 3, 42218,
                                                                       42239, 58248, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 79436, 3, 42239,
                                                                       42260, 58276, ncols,
                                                                       gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 79472, 0, 3,
                                                                       78968, 57856, 79004,
                                                                       42302, 42365, 58472,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 79580, 0, 3,
                                                                       79004, 57884, 79040,
                                                                       42365, 42428, 58556,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 79688, 0, 3,
                                                                       79040, 57912, 79076,
                                                                       42428, 42491, 58640,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 79796, 0, 3,
                                                                       79076, 57940, 79112,
                                                                       42491, 42554, 58724,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 79904, 0, 3,
                                                                       79112, 57968, 79148,
                                                                       42554, 42617, 58808,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 80012, 0, 3,
                                                                       79148, 57996, 79184,
                                                                       42617, 42680, 58892,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 80120, 0, 3,
                                                                       79220, 58108, 79256,
                                                                       42806, 42869, 59144,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 80228, 0, 3,
                                                                       79256, 58136, 79292,
                                                                       42869, 42932, 59228,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 80336, 0, 3,
                                                                       79292, 58164, 79328,
                                                                       42932, 42995, 59312,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 80444, 0, 3,
                                                                       79328, 58192, 79364,
                                                                       42995, 43058, 59396,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 80552, 0, 3,
                                                                       79364, 58220, 79400,
                                                                       43058, 43121, 59480,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 80660, 0, 3,
                                                                       79400, 58248, 79436,
                                                                       43121, 43184, 59564,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 80768, 0, 3,
                                                                       79472, 58472, 79580,
                                                                       43310, 43436, 59984,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 80984, 0, 3,
                                                                       79580, 58556, 79688,
                                                                       43436, 43562, 60152,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 81200, 0, 3,
                                                                       79688, 58640, 79796,
                                                                       43562, 43688, 60320,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 81416, 0, 3,
                                                                       79796, 58724, 79904,
                                                                       43688, 43814, 60488,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 81632, 0, 3,
                                                                       79904, 58808, 80012,
                                                                       43814, 43940, 60656,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 81848, 0, 3,
                                                                       80120, 59144, 80228,
                                                                       44192, 44318, 61160,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 82064, 0, 3,
                                                                       80228, 59228, 80336,
                                                                       44318, 44444, 61328,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 82280, 0, 3,
                                                                       80336, 59312, 80444,
                                                                       44444, 44570, 61496,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 82496, 0, 3,
                                                                       80444, 59396, 80552,
                                                                       44570, 44696, 61664,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 82712, 0, 3,
                                                                       80552, 59480, 80660,
                                                                       44696, 44822, 61832,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 82928, 0, 3,
                                                                       80768, 59984, 80984,
                                                                       45074, 45284, 62560,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 83288, 0, 3,
                                                                       80984, 60152, 81200,
                                                                       45284, 45494, 62840,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 83648, 0, 3,
                                                                       81200, 60320, 81416,
                                                                       45494, 45704, 63120,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 84008, 0, 3,
                                                                       81416, 60488, 81632,
                                                                       45704, 45914, 63400,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 84368, 0, 3,
                                                                       81848, 61160, 82064,
                                                                       46334, 46544, 64240,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 84728, 0, 3,
                                                                       82064, 61328, 82280,
                                                                       46544, 46754, 64520,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 85088, 0, 3,
                                                                       82280, 61496, 82496,
                                                                       46754, 46964, 64800,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 85448, 0, 3,
                                                                       82496, 61664, 82712,
                                                                       46964, 47174, 65080,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 85808, 0, 3,
                                                                       82928, 62560, 83288,
                                                                       47594, 47909, 66200,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 86348, 0, 3,
                                                                       83288, 62840, 83648,
                                                                       47909, 48224, 66620,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 86888, 0, 3,
                                                                       83648, 63120, 84008,
                                                                       48224, 48539, 67040,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 87428, 0, 3,
                                                                       84368, 64240, 84728,
                                                                       49169, 49484, 68300,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 87968, 0, 3,
                                                                       84728, 64520, 85088,
                                                                       49484, 49799, 68720,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 88508, 0, 3,
                                                                       85088, 64800, 85448,
                                                                       49799, 50114, 69140,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 89048, 0, 3,
                                                                       85808, 66200, 86348,
                                                                       50744, 51185, 70736,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 89804, 0, 3,
                                                                       86348, 66620, 86888,
                                                                       51185, 51626, 71324,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 90560, 0, 3,
                                                                       87428, 68300, 87968,
                                                                       52508, 52949, 73088,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 91316, 0, 3,
                                                                       87968, 68720, 88508,
                                                                       52949, 53390, 73676,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 92072, 0, 3,
                                                                       89048, 70736, 89804,
                                                                       54272, 54860, 75832,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 93080, 0, 3,
                                                                       90560, 73088, 91316,
                                                                       56036, 56624, 78184,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 94088, 3, 57800,
                                                                       57828, 78968, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 94133, 3, 57828,
                                                                       57856, 79004, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 94178, 3, 57856,
                                                                       57884, 79040, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 94223, 3, 57884,
                                                                       57912, 79076, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 94268, 3, 57912,
                                                                       57940, 79112, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 94313, 3, 57940,
                                                                       57968, 79148, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 94358, 3, 57968,
                                                                       57996, 79184, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 94403, 3, 58052,
                                                                       58080, 79220, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 94448, 3, 58080,
                                                                       58108, 79256, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 94493, 3, 58108,
                                                                       58136, 79292, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 94538, 3, 58136,
                                                                       58164, 79328, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 94583, 3, 58164,
                                                                       58192, 79364, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 94628, 3, 58192,
                                                                       58220, 79400, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 94673, 3, 58220,
                                                                       58248, 79436, ncols,
                                                                       gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 94718, 0, 3,
                                                                       94088, 78968, 94133,
                                                                       58304, 58388, 79472,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 94853, 0, 3,
                                                                       94133, 79004, 94178,
                                                                       58388, 58472, 79580,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 94988, 0, 3,
                                                                       94178, 79040, 94223,
                                                                       58472, 58556, 79688,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 95123, 0, 3,
                                                                       94223, 79076, 94268,
                                                                       58556, 58640, 79796,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 95258, 0, 3,
                                                                       94268, 79112, 94313,
                                                                       58640, 58724, 79904,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 95393, 0, 3,
                                                                       94313, 79148, 94358,
                                                                       58724, 58808, 80012,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 95528, 0, 3,
                                                                       94403, 79220, 94448,
                                                                       58976, 59060, 80120,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 95663, 0, 3,
                                                                       94448, 79256, 94493,
                                                                       59060, 59144, 80228,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 95798, 0, 3,
                                                                       94493, 79292, 94538,
                                                                       59144, 59228, 80336,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 95933, 0, 3,
                                                                       94538, 79328, 94583,
                                                                       59228, 59312, 80444,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 96068, 0, 3,
                                                                       94583, 79364, 94628,
                                                                       59312, 59396, 80552,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 96203, 0, 3,
                                                                       94628, 79400, 94673,
                                                                       59396, 59480, 80660,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 96338, 0, 3,
                                                                       94718, 79472, 94853,
                                                                       59648, 59816, 80768,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 96608, 0, 3,
                                                                       94853, 79580, 94988,
                                                                       59816, 59984, 80984,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 96878, 0, 3,
                                                                       94988, 79688, 95123,
                                                                       59984, 60152, 81200,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 97148, 0, 3,
                                                                       95123, 79796, 95258,
                                                                       60152, 60320, 81416,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 97418, 0, 3,
                                                                       95258, 79904, 95393,
                                                                       60320, 60488, 81632,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 97688, 0, 3,
                                                                       95528, 80120, 95663,
                                                                       60824, 60992, 81848,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 97958, 0, 3,
                                                                       95663, 80228, 95798,
                                                                       60992, 61160, 82064,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 98228, 0, 3,
                                                                       95798, 80336, 95933,
                                                                       61160, 61328, 82280,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 98498, 0, 3,
                                                                       95933, 80444, 96068,
                                                                       61328, 61496, 82496,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 98768, 0, 3,
                                                                       96068, 80552, 96203,
                                                                       61496, 61664, 82712,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 99038, 0, 3,
                                                                       96338, 80768, 96608,
                                                                       62000, 62280, 82928,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 99488, 0, 3,
                                                                       96608, 80984, 96878,
                                                                       62280, 62560, 83288,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 99938, 0, 3,
                                                                       96878, 81200, 97148,
                                                                       62560, 62840, 83648,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 100388, 0, 3,
                                                                       97148, 81416, 97418,
                                                                       62840, 63120, 84008,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 100838, 0, 3,
                                                                       97688, 81848, 97958,
                                                                       63680, 63960, 84368,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 101288, 0, 3,
                                                                       97958, 82064, 98228,
                                                                       63960, 64240, 84728,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 101738, 0, 3,
                                                                       98228, 82280, 98498,
                                                                       64240, 64520, 85088,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 102188, 0, 3,
                                                                       98498, 82496, 98768,
                                                                       64520, 64800, 85448,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 102638, 0, 3,
                                                                       99038, 82928, 99488,
                                                                       65360, 65780, 85808,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 103313, 0, 3,
                                                                       99488, 83288, 99938,
                                                                       65780, 66200, 86348,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 103988, 0, 3,
                                                                       99938, 83648, 100388,
                                                                       66200, 66620, 86888,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 104663, 0, 3,
                                                                       100838, 84368, 101288,
                                                                       67460, 67880, 87428,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 105338, 0, 3,
                                                                       101288, 84728, 101738,
                                                                       67880, 68300, 87968,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 106013, 0, 3,
                                                                       101738, 85088, 102188,
                                                                       68300, 68720, 88508,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 106688, 0, 3,
                                                                       102638, 85808, 103313,
                                                                       69560, 70148, 89048,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 107633, 0, 3,
                                                                       103313, 86348, 103988,
                                                                       70148, 70736, 89804,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 108578, 0, 3,
                                                                       104663, 87428, 105338,
                                                                       71912, 72500, 90560,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 109523, 0, 3,
                                                                       105338, 87968, 106013,
                                                                       72500, 73088, 91316,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 110468, 0, 3,
                                                                       106688, 89048, 107633,
                                                                       74264, 75048, 92072,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 111728, 0, 3,
                                                                       108578, 90560, 109523,
                                                                       76616, 77400, 93080,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 112988, 102638, 675, ncols);

                    simdfunc::contract_primitives(buffer, 113918, 104663, 675, ncols);

                    simdfunc::contract_primitives(buffer, 114848, 106688, 945, ncols);

                    simdfunc::contract_primitives(buffer, 116150, 108578, 945, ncols);

                    simdfunc::contract_primitives(buffer, 117452, 110468, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 119188, 111728, 1260, ncols);
                }
            }
        }

        simdtrf::transform_l_inner(buffer, 113663, 112988, 15, 1, nmax);

        simdtrf::transform_l_inner(buffer, 114593, 113918, 15, 1, nmax);

        simdtrf::transform_l_inner(buffer, 115793, 114848, 21, 1, nmax);

        simdtrf::transform_l_inner(buffer, 117095, 116150, 21, 1, nmax);

        simdtrf::transform_l_inner(buffer, 118712, 117452, 28, 1, nmax);

        simdtrf::transform_l_inner(buffer, 120448, 119188, 28, 1, nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 120924, 113663, 115793, 17, nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 121689, 114593, 117095, 17, nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 122454, 115793, 118712, 17, nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 123525, 117095, 120448, 17, nmax);

        simdtrf::compute_hrr_dg(buffer, coordinates, 124596, 120924, 122454, 17, nmax);

        simdtrf::compute_hrr_dg(buffer, coordinates, 126126, 121689, 123525, 17, nmax);

        simdtrf::transform_g_inner(buffer, 127656, 126126, 6, 17, nmax);

        simdtrf::transform_d_outer(values + n * npairs, nvalues, buffer, 127656, 153, nmax);

        simdtrf::transform_g_inner(buffer, 127656, 124596, 6, 17, nmax);

        simdtrf::transform_d_outer(values + 765 * nvalues + n * npairs, nvalues, buffer, 127656,
                                   153, nmax);
    }

    for (size_t m = 0; m < 1530; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
