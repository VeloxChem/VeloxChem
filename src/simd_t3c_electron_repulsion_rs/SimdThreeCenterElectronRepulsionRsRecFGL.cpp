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


#include "SimdThreeCenterElectronRepulsionRsRecFGL.hpp"

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
#include "SimdTransferDG.hpp"
#include "SimdTransferDH.hpp"
#include "SimdTransferFG.hpp"
#include "SimdTransferPG.hpp"
#include "SimdTransferPH.hpp"
#include "SimdTransferPI.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformL.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_fgl_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_fgl_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 204210, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 2142 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 204210, 171308, 11788, dimensions);

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

                    simdfunc::compute_full_t3c_erf_boys_function(buffer, coordinates, 6, 3, 15,
                                                                 ncols, fj, i * nprim_b + j, fq,
                                                                 omega);

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 23, 3, 15,
                                                             ncols, fj, i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 40, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 43, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 46, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 49, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 52, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 55, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 58, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 61, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 64, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 67, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 70, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 73, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 76, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 79, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 82, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 85, 0, 3, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 88, 0, 3, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 91, 0, 3, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 94, 0, 3, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 97, 0, 3, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 100, 0, 3, 29, 30,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 103, 0, 3, 30, 31,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 106, 0, 3, 31, 32,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 109, 0, 3, 32, 33,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 112, 0, 3, 33, 34,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 115, 0, 3, 34, 35,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 118, 0, 3, 35, 36,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 121, 0, 3, 36, 37,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 124, 0, 3, 37, 38,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 127, 0, 3, 38, 39,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 130, 0, 3, 7, 8,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 136, 0, 3, 8, 9,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 142, 0, 3, 9, 10,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 148, 0, 3, 10, 11,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 154, 0, 3, 11, 12,
                                                                       52, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 160, 0, 3, 12, 13,
                                                                       55, 58, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 166, 0, 3, 13, 14,
                                                                       58, 61, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 172, 0, 3, 14, 15,
                                                                       61, 64, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 178, 0, 3, 15, 16,
                                                                       64, 67, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 184, 0, 3, 16, 17,
                                                                       67, 70, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 190, 0, 3, 17, 18,
                                                                       70, 73, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 196, 0, 3, 18, 19,
                                                                       73, 76, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 202, 0, 3, 19, 20,
                                                                       76, 79, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 208, 0, 3, 20, 21,
                                                                       79, 82, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 214, 0, 3, 24, 25,
                                                                       85, 88, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 220, 0, 3, 25, 26,
                                                                       88, 91, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 226, 0, 3, 26, 27,
                                                                       91, 94, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 232, 0, 3, 27, 28,
                                                                       94, 97, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 238, 0, 3, 28, 29,
                                                                       97, 100, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 244, 0, 3, 29, 30,
                                                                       100, 103, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 250, 0, 3, 30, 31,
                                                                       103, 106, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 256, 0, 3, 31, 32,
                                                                       106, 109, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 262, 0, 3, 32, 33,
                                                                       109, 112, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 268, 0, 3, 33, 34,
                                                                       112, 115, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 274, 0, 3, 34, 35,
                                                                       115, 118, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 280, 0, 3, 35, 36,
                                                                       118, 121, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 286, 0, 3, 36, 37,
                                                                       121, 124, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 292, 0, 3, 37, 38,
                                                                       124, 127, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 298, 0, 3, 40, 43,
                                                                       130, 136, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 308, 0, 3, 43, 46,
                                                                       136, 142, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 318, 0, 3, 46, 49,
                                                                       142, 148, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 328, 0, 3, 49, 52,
                                                                       148, 154, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 338, 0, 3, 52, 55,
                                                                       154, 160, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 348, 0, 3, 55, 58,
                                                                       160, 166, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 358, 0, 3, 58, 61,
                                                                       166, 172, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 368, 0, 3, 61, 64,
                                                                       172, 178, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 378, 0, 3, 64, 67,
                                                                       178, 184, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 388, 0, 3, 67, 70,
                                                                       184, 190, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 398, 0, 3, 70, 73,
                                                                       190, 196, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 408, 0, 3, 73, 76,
                                                                       196, 202, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 418, 0, 3, 76, 79,
                                                                       202, 208, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 428, 0, 3, 85, 88,
                                                                       214, 220, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 438, 0, 3, 88, 91,
                                                                       220, 226, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 448, 0, 3, 91, 94,
                                                                       226, 232, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 458, 0, 3, 94, 97,
                                                                       232, 238, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 468, 0, 3, 97,
                                                                       100, 238, 244, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 478, 0, 3, 100,
                                                                       103, 244, 250, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 488, 0, 3, 103,
                                                                       106, 250, 256, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 498, 0, 3, 106,
                                                                       109, 256, 262, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 508, 0, 3, 109,
                                                                       112, 262, 268, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 518, 0, 3, 112,
                                                                       115, 268, 274, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 528, 0, 3, 115,
                                                                       118, 274, 280, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 538, 0, 3, 118,
                                                                       121, 280, 286, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 548, 0, 3, 121,
                                                                       124, 286, 292, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 558, 0, 3, 130,
                                                                       136, 298, 308, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 573, 0, 3, 136,
                                                                       142, 308, 318, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 588, 0, 3, 142,
                                                                       148, 318, 328, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 603, 0, 3, 148,
                                                                       154, 328, 338, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 618, 0, 3, 154,
                                                                       160, 338, 348, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 633, 0, 3, 160,
                                                                       166, 348, 358, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 648, 0, 3, 166,
                                                                       172, 358, 368, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 663, 0, 3, 172,
                                                                       178, 368, 378, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 678, 0, 3, 178,
                                                                       184, 378, 388, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 693, 0, 3, 184,
                                                                       190, 388, 398, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 708, 0, 3, 190,
                                                                       196, 398, 408, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 723, 0, 3, 196,
                                                                       202, 408, 418, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 738, 0, 3, 214,
                                                                       220, 428, 438, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 753, 0, 3, 220,
                                                                       226, 438, 448, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 768, 0, 3, 226,
                                                                       232, 448, 458, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 783, 0, 3, 232,
                                                                       238, 458, 468, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 798, 0, 3, 238,
                                                                       244, 468, 478, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 813, 0, 3, 244,
                                                                       250, 478, 488, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 828, 0, 3, 250,
                                                                       256, 488, 498, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 843, 0, 3, 256,
                                                                       262, 498, 508, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 858, 0, 3, 262,
                                                                       268, 508, 518, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 873, 0, 3, 268,
                                                                       274, 518, 528, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 888, 0, 3, 274,
                                                                       280, 528, 538, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 903, 0, 3, 280,
                                                                       286, 538, 548, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 918, 0, 3, 298,
                                                                       308, 558, 573, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 939, 0, 3, 308,
                                                                       318, 573, 588, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 960, 0, 3, 318,
                                                                       328, 588, 603, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 981, 0, 3, 328,
                                                                       338, 603, 618, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1002, 0, 3, 338,
                                                                       348, 618, 633, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1023, 0, 3, 348,
                                                                       358, 633, 648, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1044, 0, 3, 358,
                                                                       368, 648, 663, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1065, 0, 3, 368,
                                                                       378, 663, 678, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1086, 0, 3, 378,
                                                                       388, 678, 693, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1107, 0, 3, 388,
                                                                       398, 693, 708, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1128, 0, 3, 398,
                                                                       408, 708, 723, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1149, 0, 3, 428,
                                                                       438, 738, 753, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1170, 0, 3, 438,
                                                                       448, 753, 768, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1191, 0, 3, 448,
                                                                       458, 768, 783, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1212, 0, 3, 458,
                                                                       468, 783, 798, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1233, 0, 3, 468,
                                                                       478, 798, 813, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1254, 0, 3, 478,
                                                                       488, 813, 828, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1275, 0, 3, 488,
                                                                       498, 828, 843, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1296, 0, 3, 498,
                                                                       508, 843, 858, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1317, 0, 3, 508,
                                                                       518, 858, 873, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1338, 0, 3, 518,
                                                                       528, 873, 888, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1359, 0, 3, 528,
                                                                       538, 888, 903, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1380, 0, 3, 558,
                                                                       573, 918, 939, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1408, 0, 3, 573,
                                                                       588, 939, 960, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1436, 0, 3, 588,
                                                                       603, 960, 981, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1464, 0, 3, 603,
                                                                       618, 981, 1002, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1492, 0, 3, 618,
                                                                       633, 1002, 1023, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1520, 0, 3, 633,
                                                                       648, 1023, 1044, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1548, 0, 3, 648,
                                                                       663, 1044, 1065, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1576, 0, 3, 663,
                                                                       678, 1065, 1086, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1604, 0, 3, 678,
                                                                       693, 1086, 1107, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1632, 0, 3, 693,
                                                                       708, 1107, 1128, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1660, 0, 3, 738,
                                                                       753, 1149, 1170, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1688, 0, 3, 753,
                                                                       768, 1170, 1191, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1716, 0, 3, 768,
                                                                       783, 1191, 1212, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1744, 0, 3, 783,
                                                                       798, 1212, 1233, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1772, 0, 3, 798,
                                                                       813, 1233, 1254, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1800, 0, 3, 813,
                                                                       828, 1254, 1275, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1828, 0, 3, 828,
                                                                       843, 1275, 1296, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1856, 0, 3, 843,
                                                                       858, 1296, 1317, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1884, 0, 3, 858,
                                                                       873, 1317, 1338, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1912, 0, 3, 873,
                                                                       888, 1338, 1359, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1940, 0, 3, 918,
                                                                       939, 1380, 1408, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1976, 0, 3, 939,
                                                                       960, 1408, 1436, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2012, 0, 3, 960,
                                                                       981, 1436, 1464, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2048, 0, 3, 981,
                                                                       1002, 1464, 1492, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2084, 0, 3, 1002,
                                                                       1023, 1492, 1520, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2120, 0, 3, 1023,
                                                                       1044, 1520, 1548, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2156, 0, 3, 1044,
                                                                       1065, 1548, 1576, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2192, 0, 3, 1065,
                                                                       1086, 1576, 1604, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2228, 0, 3, 1086,
                                                                       1107, 1604, 1632, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2264, 0, 3, 1149,
                                                                       1170, 1660, 1688, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2300, 0, 3, 1170,
                                                                       1191, 1688, 1716, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2336, 0, 3, 1191,
                                                                       1212, 1716, 1744, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2372, 0, 3, 1212,
                                                                       1233, 1744, 1772, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2408, 0, 3, 1233,
                                                                       1254, 1772, 1800, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2444, 0, 3, 1254,
                                                                       1275, 1800, 1828, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2480, 0, 3, 1275,
                                                                       1296, 1828, 1856, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2516, 0, 3, 1296,
                                                                       1317, 1856, 1884, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2552, 0, 3, 1317,
                                                                       1338, 1884, 1912, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2588, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2591, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2594, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2597, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2600, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2603, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2606, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2609, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2612, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2615, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2618, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2621, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2624, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2627, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2630, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2633, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2636, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2639, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2642, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2645, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2648, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2651, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2654, 3, 34,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2657, 3, 35,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2660, 3, 36,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2663, 3, 37,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2666, 3, 38,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2669, 3, 39,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2672, 3, 9, 46,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2681, 3, 10, 49,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2690, 3, 11, 52,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2699, 3, 12, 55,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2708, 3, 13, 58,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2717, 3, 14, 61,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2726, 3, 15, 64,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2735, 3, 16, 67,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2744, 3, 17, 70,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2753, 3, 18, 73,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2762, 3, 19, 76,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2771, 3, 20, 79,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2780, 3, 21, 82,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2789, 3, 26, 91,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2798, 3, 27, 94,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2807, 3, 28, 97,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2816, 3, 29, 100,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2825, 3, 30, 103,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2834, 3, 31, 106,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2843, 3, 32, 109,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2852, 3, 33, 112,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2861, 3, 34, 115,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2870, 3, 35, 118,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2879, 3, 36, 121,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2888, 3, 37, 124,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2897, 3, 38, 127,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2906, 3, 46, 142,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2924, 3, 49, 148,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2942, 3, 52, 154,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2960, 3, 55, 160,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2978, 3, 58, 166,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2996, 3, 61, 172,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3014, 3, 64, 178,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3032, 3, 67, 184,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3050, 3, 70, 190,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3068, 3, 73, 196,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3086, 3, 76, 202,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3104, 3, 79, 208,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3122, 3, 91, 226,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3140, 3, 94, 232,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3158, 3, 97, 238,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3176, 3, 100, 244,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3194, 3, 103, 250,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3212, 3, 106, 256,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3230, 3, 109, 262,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3248, 3, 112, 268,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3266, 3, 115, 274,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3284, 3, 118, 280,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3302, 3, 121, 286,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3320, 3, 124, 292,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3338, 3, 142, 318,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3368, 3, 148, 328,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3398, 3, 154, 338,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3428, 3, 160, 348,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3458, 3, 166, 358,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3488, 3, 172, 368,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3518, 3, 178, 378,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3548, 3, 184, 388,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3578, 3, 190, 398,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3608, 3, 196, 408,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3638, 3, 202, 418,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3668, 3, 226, 448,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3698, 3, 232, 458,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3728, 3, 238, 468,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3758, 3, 244, 478,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3788, 3, 250, 488,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3818, 3, 256, 498,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3848, 3, 262, 508,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3878, 3, 268, 518,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3908, 3, 274, 528,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3938, 3, 280, 538,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3968, 3, 286, 548,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3998, 3, 318, 588,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4043, 3, 328, 603,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4088, 3, 338, 618,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4133, 3, 348, 633,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4178, 3, 358, 648,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4223, 3, 368, 663,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4268, 3, 378, 678,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4313, 3, 388, 693,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4358, 3, 398, 708,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4403, 3, 408, 723,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4448, 3, 448, 768,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4493, 3, 458, 783,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4538, 3, 468, 798,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4583, 3, 478, 813,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4628, 3, 488, 828,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4673, 3, 498, 843,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4718, 3, 508, 858,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4763, 3, 518, 873,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4808, 3, 528, 888,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4853, 3, 538, 903,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4898, 3, 588, 960,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4961, 3, 603, 981,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5024, 3, 618,
                                                                       1002, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5087, 3, 633,
                                                                       1023, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5150, 3, 648,
                                                                       1044, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5213, 3, 663,
                                                                       1065, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5276, 3, 678,
                                                                       1086, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5339, 3, 693,
                                                                       1107, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5402, 3, 708,
                                                                       1128, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5465, 3, 768,
                                                                       1191, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5528, 3, 783,
                                                                       1212, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5591, 3, 798,
                                                                       1233, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5654, 3, 813,
                                                                       1254, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5717, 3, 828,
                                                                       1275, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5780, 3, 843,
                                                                       1296, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5843, 3, 858,
                                                                       1317, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5906, 3, 873,
                                                                       1338, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5969, 3, 888,
                                                                       1359, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6032, 3, 960,
                                                                       1436, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6116, 3, 981,
                                                                       1464, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6200, 3, 1002,
                                                                       1492, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6284, 3, 1023,
                                                                       1520, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6368, 3, 1044,
                                                                       1548, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6452, 3, 1065,
                                                                       1576, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6536, 3, 1086,
                                                                       1604, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6620, 3, 1107,
                                                                       1632, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6704, 3, 1191,
                                                                       1716, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6788, 3, 1212,
                                                                       1744, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6872, 3, 1233,
                                                                       1772, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6956, 3, 1254,
                                                                       1800, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7040, 3, 1275,
                                                                       1828, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7124, 3, 1296,
                                                                       1856, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7208, 3, 1317,
                                                                       1884, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7292, 3, 1338,
                                                                       1912, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7376, 3, 1436,
                                                                       2012, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7484, 3, 1464,
                                                                       2048, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7592, 3, 1492,
                                                                       2084, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7700, 3, 1520,
                                                                       2120, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7808, 3, 1548,
                                                                       2156, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7916, 3, 1576,
                                                                       2192, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8024, 3, 1604,
                                                                       2228, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8132, 3, 1716,
                                                                       2336, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8240, 3, 1744,
                                                                       2372, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8348, 3, 1772,
                                                                       2408, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8456, 3, 1800,
                                                                       2444, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8564, 3, 1828,
                                                                       2480, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8672, 3, 1856,
                                                                       2516, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8780, 3, 1884,
                                                                       2552, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8888, 3, 7, 8,
                                                                       2588, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8894, 3, 8, 9,
                                                                       2591, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8900, 3, 9, 10,
                                                                       2594, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8906, 3, 10, 11,
                                                                       2597, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8912, 3, 11, 12,
                                                                       2600, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8918, 3, 12, 13,
                                                                       2603, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8924, 3, 13, 14,
                                                                       2606, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8930, 3, 14, 15,
                                                                       2609, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8936, 3, 15, 16,
                                                                       2612, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8942, 3, 16, 17,
                                                                       2615, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8948, 3, 17, 18,
                                                                       2618, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8954, 3, 18, 19,
                                                                       2621, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8960, 3, 19, 20,
                                                                       2624, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8966, 3, 20, 21,
                                                                       2627, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8972, 3, 24, 25,
                                                                       2630, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8978, 3, 25, 26,
                                                                       2633, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8984, 3, 26, 27,
                                                                       2636, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8990, 3, 27, 28,
                                                                       2639, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8996, 3, 28, 29,
                                                                       2642, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9002, 3, 29, 30,
                                                                       2645, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9008, 3, 30, 31,
                                                                       2648, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9014, 3, 31, 32,
                                                                       2651, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9020, 3, 32, 33,
                                                                       2654, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9026, 3, 33, 34,
                                                                       2657, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9032, 3, 34, 35,
                                                                       2660, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9038, 3, 35, 36,
                                                                       2663, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9044, 3, 36, 37,
                                                                       2666, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9050, 3, 37, 38,
                                                                       2669, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9056, 0, 3, 8888,
                                                                       2588, 8894, 2672, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9074, 0, 3, 8894,
                                                                       2591, 8900, 2681, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9092, 0, 3, 8900,
                                                                       2594, 8906, 2690, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9110, 0, 3, 8906,
                                                                       2597, 8912, 2699, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9128, 0, 3, 8912,
                                                                       2600, 8918, 2708, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9146, 0, 3, 8918,
                                                                       2603, 8924, 2717, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9164, 0, 3, 8924,
                                                                       2606, 8930, 2726, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9182, 0, 3, 8930,
                                                                       2609, 8936, 2735, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9200, 0, 3, 8936,
                                                                       2612, 8942, 2744, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9218, 0, 3, 8942,
                                                                       2615, 8948, 2753, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9236, 0, 3, 8948,
                                                                       2618, 8954, 2762, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9254, 0, 3, 8954,
                                                                       2621, 8960, 2771, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9272, 0, 3, 8960,
                                                                       2624, 8966, 2780, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9290, 0, 3, 8972,
                                                                       2630, 8978, 2789, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9308, 0, 3, 8978,
                                                                       2633, 8984, 2798, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9326, 0, 3, 8984,
                                                                       2636, 8990, 2807, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9344, 0, 3, 8990,
                                                                       2639, 8996, 2816, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9362, 0, 3, 8996,
                                                                       2642, 9002, 2825, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9380, 0, 3, 9002,
                                                                       2645, 9008, 2834, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9398, 0, 3, 9008,
                                                                       2648, 9014, 2843, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9416, 0, 3, 9014,
                                                                       2651, 9020, 2852, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9434, 0, 3, 9020,
                                                                       2654, 9026, 2861, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9452, 0, 3, 9026,
                                                                       2657, 9032, 2870, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9470, 0, 3, 9032,
                                                                       2660, 9038, 2879, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9488, 0, 3, 9038,
                                                                       2663, 9044, 2888, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9506, 0, 3, 9044,
                                                                       2666, 9050, 2897, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9524, 0, 3, 9056,
                                                                       2672, 9074, 130, 136,
                                                                       2906, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9560, 0, 3, 9074,
                                                                       2681, 9092, 136, 142,
                                                                       2924, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9596, 0, 3, 9092,
                                                                       2690, 9110, 142, 148,
                                                                       2942, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9632, 0, 3, 9110,
                                                                       2699, 9128, 148, 154,
                                                                       2960, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9668, 0, 3, 9128,
                                                                       2708, 9146, 154, 160,
                                                                       2978, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9704, 0, 3, 9146,
                                                                       2717, 9164, 160, 166,
                                                                       2996, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9740, 0, 3, 9164,
                                                                       2726, 9182, 166, 172,
                                                                       3014, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9776, 0, 3, 9182,
                                                                       2735, 9200, 172, 178,
                                                                       3032, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9812, 0, 3, 9200,
                                                                       2744, 9218, 178, 184,
                                                                       3050, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9848, 0, 3, 9218,
                                                                       2753, 9236, 184, 190,
                                                                       3068, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9884, 0, 3, 9236,
                                                                       2762, 9254, 190, 196,
                                                                       3086, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9920, 0, 3, 9254,
                                                                       2771, 9272, 196, 202,
                                                                       3104, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9956, 0, 3, 9290,
                                                                       2789, 9308, 214, 220,
                                                                       3122, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9992, 0, 3, 9308,
                                                                       2798, 9326, 220, 226,
                                                                       3140, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 10028, 0, 3, 9326,
                                                                       2807, 9344, 226, 232,
                                                                       3158, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 10064, 0, 3, 9344,
                                                                       2816, 9362, 232, 238,
                                                                       3176, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 10100, 0, 3, 9362,
                                                                       2825, 9380, 238, 244,
                                                                       3194, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 10136, 0, 3, 9380,
                                                                       2834, 9398, 244, 250,
                                                                       3212, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 10172, 0, 3, 9398,
                                                                       2843, 9416, 250, 256,
                                                                       3230, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 10208, 0, 3, 9416,
                                                                       2852, 9434, 256, 262,
                                                                       3248, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 10244, 0, 3, 9434,
                                                                       2861, 9452, 262, 268,
                                                                       3266, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 10280, 0, 3, 9452,
                                                                       2870, 9470, 268, 274,
                                                                       3284, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 10316, 0, 3, 9470,
                                                                       2879, 9488, 274, 280,
                                                                       3302, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 10352, 0, 3, 9488,
                                                                       2888, 9506, 280, 286,
                                                                       3320, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10388, 0, 3, 9524,
                                                                       2906, 9560, 298, 308,
                                                                       3338, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10448, 0, 3, 9560,
                                                                       2924, 9596, 308, 318,
                                                                       3368, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10508, 0, 3, 9596,
                                                                       2942, 9632, 318, 328,
                                                                       3398, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10568, 0, 3, 9632,
                                                                       2960, 9668, 328, 338,
                                                                       3428, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10628, 0, 3, 9668,
                                                                       2978, 9704, 338, 348,
                                                                       3458, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10688, 0, 3, 9704,
                                                                       2996, 9740, 348, 358,
                                                                       3488, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10748, 0, 3, 9740,
                                                                       3014, 9776, 358, 368,
                                                                       3518, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10808, 0, 3, 9776,
                                                                       3032, 9812, 368, 378,
                                                                       3548, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10868, 0, 3, 9812,
                                                                       3050, 9848, 378, 388,
                                                                       3578, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10928, 0, 3, 9848,
                                                                       3068, 9884, 388, 398,
                                                                       3608, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10988, 0, 3, 9884,
                                                                       3086, 9920, 398, 408,
                                                                       3638, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 11048, 0, 3, 9956,
                                                                       3122, 9992, 428, 438,
                                                                       3668, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 11108, 0, 3, 9992,
                                                                       3140, 10028, 438, 448,
                                                                       3698, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 11168, 0, 3,
                                                                       10028, 3158, 10064, 448,
                                                                       458, 3728, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 11228, 0, 3,
                                                                       10064, 3176, 10100, 458,
                                                                       468, 3758, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 11288, 0, 3,
                                                                       10100, 3194, 10136, 468,
                                                                       478, 3788, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 11348, 0, 3,
                                                                       10136, 3212, 10172, 478,
                                                                       488, 3818, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 11408, 0, 3,
                                                                       10172, 3230, 10208, 488,
                                                                       498, 3848, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 11468, 0, 3,
                                                                       10208, 3248, 10244, 498,
                                                                       508, 3878, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 11528, 0, 3,
                                                                       10244, 3266, 10280, 508,
                                                                       518, 3908, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 11588, 0, 3,
                                                                       10280, 3284, 10316, 518,
                                                                       528, 3938, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 11648, 0, 3,
                                                                       10316, 3302, 10352, 528,
                                                                       538, 3968, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11708, 0, 3,
                                                                       10388, 3338, 10448, 558,
                                                                       573, 3998, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11798, 0, 3,
                                                                       10448, 3368, 10508, 573,
                                                                       588, 4043, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11888, 0, 3,
                                                                       10508, 3398, 10568, 588,
                                                                       603, 4088, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11978, 0, 3,
                                                                       10568, 3428, 10628, 603,
                                                                       618, 4133, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12068, 0, 3,
                                                                       10628, 3458, 10688, 618,
                                                                       633, 4178, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12158, 0, 3,
                                                                       10688, 3488, 10748, 633,
                                                                       648, 4223, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12248, 0, 3,
                                                                       10748, 3518, 10808, 648,
                                                                       663, 4268, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12338, 0, 3,
                                                                       10808, 3548, 10868, 663,
                                                                       678, 4313, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12428, 0, 3,
                                                                       10868, 3578, 10928, 678,
                                                                       693, 4358, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12518, 0, 3,
                                                                       10928, 3608, 10988, 693,
                                                                       708, 4403, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12608, 0, 3,
                                                                       11048, 3668, 11108, 738,
                                                                       753, 4448, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12698, 0, 3,
                                                                       11108, 3698, 11168, 753,
                                                                       768, 4493, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12788, 0, 3,
                                                                       11168, 3728, 11228, 768,
                                                                       783, 4538, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12878, 0, 3,
                                                                       11228, 3758, 11288, 783,
                                                                       798, 4583, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12968, 0, 3,
                                                                       11288, 3788, 11348, 798,
                                                                       813, 4628, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 13058, 0, 3,
                                                                       11348, 3818, 11408, 813,
                                                                       828, 4673, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 13148, 0, 3,
                                                                       11408, 3848, 11468, 828,
                                                                       843, 4718, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 13238, 0, 3,
                                                                       11468, 3878, 11528, 843,
                                                                       858, 4763, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 13328, 0, 3,
                                                                       11528, 3908, 11588, 858,
                                                                       873, 4808, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 13418, 0, 3,
                                                                       11588, 3938, 11648, 873,
                                                                       888, 4853, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 13508, 0, 3,
                                                                       11708, 3998, 11798, 918,
                                                                       939, 4898, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 13634, 0, 3,
                                                                       11798, 4043, 11888, 939,
                                                                       960, 4961, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 13760, 0, 3,
                                                                       11888, 4088, 11978, 960,
                                                                       981, 5024, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 13886, 0, 3,
                                                                       11978, 4133, 12068, 981,
                                                                       1002, 5087, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 14012, 0, 3,
                                                                       12068, 4178, 12158, 1002,
                                                                       1023, 5150, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 14138, 0, 3,
                                                                       12158, 4223, 12248, 1023,
                                                                       1044, 5213, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 14264, 0, 3,
                                                                       12248, 4268, 12338, 1044,
                                                                       1065, 5276, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 14390, 0, 3,
                                                                       12338, 4313, 12428, 1065,
                                                                       1086, 5339, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 14516, 0, 3,
                                                                       12428, 4358, 12518, 1086,
                                                                       1107, 5402, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 14642, 0, 3,
                                                                       12608, 4448, 12698, 1149,
                                                                       1170, 5465, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 14768, 0, 3,
                                                                       12698, 4493, 12788, 1170,
                                                                       1191, 5528, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 14894, 0, 3,
                                                                       12788, 4538, 12878, 1191,
                                                                       1212, 5591, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 15020, 0, 3,
                                                                       12878, 4583, 12968, 1212,
                                                                       1233, 5654, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 15146, 0, 3,
                                                                       12968, 4628, 13058, 1233,
                                                                       1254, 5717, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 15272, 0, 3,
                                                                       13058, 4673, 13148, 1254,
                                                                       1275, 5780, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 15398, 0, 3,
                                                                       13148, 4718, 13238, 1275,
                                                                       1296, 5843, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 15524, 0, 3,
                                                                       13238, 4763, 13328, 1296,
                                                                       1317, 5906, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 15650, 0, 3,
                                                                       13328, 4808, 13418, 1317,
                                                                       1338, 5969, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 15776, 0, 3,
                                                                       13508, 4898, 13634, 1380,
                                                                       1408, 6032, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 15944, 0, 3,
                                                                       13634, 4961, 13760, 1408,
                                                                       1436, 6116, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 16112, 0, 3,
                                                                       13760, 5024, 13886, 1436,
                                                                       1464, 6200, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 16280, 0, 3,
                                                                       13886, 5087, 14012, 1464,
                                                                       1492, 6284, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 16448, 0, 3,
                                                                       14012, 5150, 14138, 1492,
                                                                       1520, 6368, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 16616, 0, 3,
                                                                       14138, 5213, 14264, 1520,
                                                                       1548, 6452, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 16784, 0, 3,
                                                                       14264, 5276, 14390, 1548,
                                                                       1576, 6536, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 16952, 0, 3,
                                                                       14390, 5339, 14516, 1576,
                                                                       1604, 6620, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 17120, 0, 3,
                                                                       14642, 5465, 14768, 1660,
                                                                       1688, 6704, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 17288, 0, 3,
                                                                       14768, 5528, 14894, 1688,
                                                                       1716, 6788, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 17456, 0, 3,
                                                                       14894, 5591, 15020, 1716,
                                                                       1744, 6872, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 17624, 0, 3,
                                                                       15020, 5654, 15146, 1744,
                                                                       1772, 6956, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 17792, 0, 3,
                                                                       15146, 5717, 15272, 1772,
                                                                       1800, 7040, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 17960, 0, 3,
                                                                       15272, 5780, 15398, 1800,
                                                                       1828, 7124, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 18128, 0, 3,
                                                                       15398, 5843, 15524, 1828,
                                                                       1856, 7208, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 18296, 0, 3,
                                                                       15524, 5906, 15650, 1856,
                                                                       1884, 7292, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 18464, 0, 3,
                                                                       15776, 6032, 15944, 1940,
                                                                       1976, 7376, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 18680, 0, 3,
                                                                       15944, 6116, 16112, 1976,
                                                                       2012, 7484, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 18896, 0, 3,
                                                                       16112, 6200, 16280, 2012,
                                                                       2048, 7592, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 19112, 0, 3,
                                                                       16280, 6284, 16448, 2048,
                                                                       2084, 7700, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 19328, 0, 3,
                                                                       16448, 6368, 16616, 2084,
                                                                       2120, 7808, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 19544, 0, 3,
                                                                       16616, 6452, 16784, 2120,
                                                                       2156, 7916, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 19760, 0, 3,
                                                                       16784, 6536, 16952, 2156,
                                                                       2192, 8024, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 19976, 0, 3,
                                                                       17120, 6704, 17288, 2264,
                                                                       2300, 8132, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 20192, 0, 3,
                                                                       17288, 6788, 17456, 2300,
                                                                       2336, 8240, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 20408, 0, 3,
                                                                       17456, 6872, 17624, 2336,
                                                                       2372, 8348, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 20624, 0, 3,
                                                                       17624, 6956, 17792, 2372,
                                                                       2408, 8456, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 20840, 0, 3,
                                                                       17792, 7040, 17960, 2408,
                                                                       2444, 8564, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 21056, 0, 3,
                                                                       17960, 7124, 18128, 2444,
                                                                       2480, 8672, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 21272, 0, 3,
                                                                       18128, 7208, 18296, 2480,
                                                                       2516, 8780, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21488, 3, 2588,
                                                                       2591, 8900, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21498, 3, 2591,
                                                                       2594, 8906, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21508, 3, 2594,
                                                                       2597, 8912, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21518, 3, 2597,
                                                                       2600, 8918, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21528, 3, 2600,
                                                                       2603, 8924, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21538, 3, 2603,
                                                                       2606, 8930, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21548, 3, 2606,
                                                                       2609, 8936, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21558, 3, 2609,
                                                                       2612, 8942, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21568, 3, 2612,
                                                                       2615, 8948, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21578, 3, 2615,
                                                                       2618, 8954, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21588, 3, 2618,
                                                                       2621, 8960, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21598, 3, 2621,
                                                                       2624, 8966, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21608, 3, 2630,
                                                                       2633, 8984, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21618, 3, 2633,
                                                                       2636, 8990, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21628, 3, 2636,
                                                                       2639, 8996, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21638, 3, 2639,
                                                                       2642, 9002, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21648, 3, 2642,
                                                                       2645, 9008, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21658, 3, 2645,
                                                                       2648, 9014, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21668, 3, 2648,
                                                                       2651, 9020, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21678, 3, 2651,
                                                                       2654, 9026, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21688, 3, 2654,
                                                                       2657, 9032, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21698, 3, 2657,
                                                                       2660, 9038, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21708, 3, 2660,
                                                                       2663, 9044, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21718, 3, 2663,
                                                                       2666, 9050, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 21728, 0, 3,
                                                                       21488, 8900, 21498, 9092,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 21758, 0, 3,
                                                                       21498, 8906, 21508, 9110,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 21788, 0, 3,
                                                                       21508, 8912, 21518, 9128,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 21818, 0, 3,
                                                                       21518, 8918, 21528, 9146,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 21848, 0, 3,
                                                                       21528, 8924, 21538, 9164,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 21878, 0, 3,
                                                                       21538, 8930, 21548, 9182,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 21908, 0, 3,
                                                                       21548, 8936, 21558, 9200,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 21938, 0, 3,
                                                                       21558, 8942, 21568, 9218,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 21968, 0, 3,
                                                                       21568, 8948, 21578, 9236,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 21998, 0, 3,
                                                                       21578, 8954, 21588, 9254,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 22028, 0, 3,
                                                                       21588, 8960, 21598, 9272,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 22058, 0, 3,
                                                                       21608, 8984, 21618, 9326,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 22088, 0, 3,
                                                                       21618, 8990, 21628, 9344,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 22118, 0, 3,
                                                                       21628, 8996, 21638, 9362,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 22148, 0, 3,
                                                                       21638, 9002, 21648, 9380,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 22178, 0, 3,
                                                                       21648, 9008, 21658, 9398,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 22208, 0, 3,
                                                                       21658, 9014, 21668, 9416,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 22238, 0, 3,
                                                                       21668, 9020, 21678, 9434,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 22268, 0, 3,
                                                                       21678, 9026, 21688, 9452,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 22298, 0, 3,
                                                                       21688, 9032, 21698, 9470,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 22328, 0, 3,
                                                                       21698, 9038, 21708, 9488,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 22358, 0, 3,
                                                                       21708, 9044, 21718, 9506,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 22388, 0, 3,
                                                                       21728, 9092, 21758, 2906,
                                                                       2924, 9596, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 22448, 0, 3,
                                                                       21758, 9110, 21788, 2924,
                                                                       2942, 9632, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 22508, 0, 3,
                                                                       21788, 9128, 21818, 2942,
                                                                       2960, 9668, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 22568, 0, 3,
                                                                       21818, 9146, 21848, 2960,
                                                                       2978, 9704, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 22628, 0, 3,
                                                                       21848, 9164, 21878, 2978,
                                                                       2996, 9740, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 22688, 0, 3,
                                                                       21878, 9182, 21908, 2996,
                                                                       3014, 9776, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 22748, 0, 3,
                                                                       21908, 9200, 21938, 3014,
                                                                       3032, 9812, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 22808, 0, 3,
                                                                       21938, 9218, 21968, 3032,
                                                                       3050, 9848, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 22868, 0, 3,
                                                                       21968, 9236, 21998, 3050,
                                                                       3068, 9884, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 22928, 0, 3,
                                                                       21998, 9254, 22028, 3068,
                                                                       3086, 9920, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 22988, 0, 3,
                                                                       22058, 9326, 22088, 3122,
                                                                       3140, 10028, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 23048, 0, 3,
                                                                       22088, 9344, 22118, 3140,
                                                                       3158, 10064, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 23108, 0, 3,
                                                                       22118, 9362, 22148, 3158,
                                                                       3176, 10100, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 23168, 0, 3,
                                                                       22148, 9380, 22178, 3176,
                                                                       3194, 10136, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 23228, 0, 3,
                                                                       22178, 9398, 22208, 3194,
                                                                       3212, 10172, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 23288, 0, 3,
                                                                       22208, 9416, 22238, 3212,
                                                                       3230, 10208, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 23348, 0, 3,
                                                                       22238, 9434, 22268, 3230,
                                                                       3248, 10244, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 23408, 0, 3,
                                                                       22268, 9452, 22298, 3248,
                                                                       3266, 10280, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 23468, 0, 3,
                                                                       22298, 9470, 22328, 3266,
                                                                       3284, 10316, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 23528, 0, 3,
                                                                       22328, 9488, 22358, 3284,
                                                                       3302, 10352, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 23588, 0, 3,
                                                                       22388, 9596, 22448, 3338,
                                                                       3368, 10508, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 23688, 0, 3,
                                                                       22448, 9632, 22508, 3368,
                                                                       3398, 10568, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 23788, 0, 3,
                                                                       22508, 9668, 22568, 3398,
                                                                       3428, 10628, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 23888, 0, 3,
                                                                       22568, 9704, 22628, 3428,
                                                                       3458, 10688, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 23988, 0, 3,
                                                                       22628, 9740, 22688, 3458,
                                                                       3488, 10748, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 24088, 0, 3,
                                                                       22688, 9776, 22748, 3488,
                                                                       3518, 10808, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 24188, 0, 3,
                                                                       22748, 9812, 22808, 3518,
                                                                       3548, 10868, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 24288, 0, 3,
                                                                       22808, 9848, 22868, 3548,
                                                                       3578, 10928, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 24388, 0, 3,
                                                                       22868, 9884, 22928, 3578,
                                                                       3608, 10988, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 24488, 0, 3,
                                                                       22988, 10028, 23048, 3668,
                                                                       3698, 11168, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 24588, 0, 3,
                                                                       23048, 10064, 23108, 3698,
                                                                       3728, 11228, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 24688, 0, 3,
                                                                       23108, 10100, 23168, 3728,
                                                                       3758, 11288, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 24788, 0, 3,
                                                                       23168, 10136, 23228, 3758,
                                                                       3788, 11348, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 24888, 0, 3,
                                                                       23228, 10172, 23288, 3788,
                                                                       3818, 11408, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 24988, 0, 3,
                                                                       23288, 10208, 23348, 3818,
                                                                       3848, 11468, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25088, 0, 3,
                                                                       23348, 10244, 23408, 3848,
                                                                       3878, 11528, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25188, 0, 3,
                                                                       23408, 10280, 23468, 3878,
                                                                       3908, 11588, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25288, 0, 3,
                                                                       23468, 10316, 23528, 3908,
                                                                       3938, 11648, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 25388, 0, 3,
                                                                       23588, 10508, 23688, 3998,
                                                                       4043, 11888, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 25538, 0, 3,
                                                                       23688, 10568, 23788, 4043,
                                                                       4088, 11978, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 25688, 0, 3,
                                                                       23788, 10628, 23888, 4088,
                                                                       4133, 12068, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 25838, 0, 3,
                                                                       23888, 10688, 23988, 4133,
                                                                       4178, 12158, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 25988, 0, 3,
                                                                       23988, 10748, 24088, 4178,
                                                                       4223, 12248, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 26138, 0, 3,
                                                                       24088, 10808, 24188, 4223,
                                                                       4268, 12338, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 26288, 0, 3,
                                                                       24188, 10868, 24288, 4268,
                                                                       4313, 12428, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 26438, 0, 3,
                                                                       24288, 10928, 24388, 4313,
                                                                       4358, 12518, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 26588, 0, 3,
                                                                       24488, 11168, 24588, 4448,
                                                                       4493, 12788, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 26738, 0, 3,
                                                                       24588, 11228, 24688, 4493,
                                                                       4538, 12878, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 26888, 0, 3,
                                                                       24688, 11288, 24788, 4538,
                                                                       4583, 12968, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 27038, 0, 3,
                                                                       24788, 11348, 24888, 4583,
                                                                       4628, 13058, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 27188, 0, 3,
                                                                       24888, 11408, 24988, 4628,
                                                                       4673, 13148, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 27338, 0, 3,
                                                                       24988, 11468, 25088, 4673,
                                                                       4718, 13238, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 27488, 0, 3,
                                                                       25088, 11528, 25188, 4718,
                                                                       4763, 13328, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 27638, 0, 3,
                                                                       25188, 11588, 25288, 4763,
                                                                       4808, 13418, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 27788, 0, 3,
                                                                       25388, 11888, 25538, 4898,
                                                                       4961, 13760, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 27998, 0, 3,
                                                                       25538, 11978, 25688, 4961,
                                                                       5024, 13886, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 28208, 0, 3,
                                                                       25688, 12068, 25838, 5024,
                                                                       5087, 14012, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 28418, 0, 3,
                                                                       25838, 12158, 25988, 5087,
                                                                       5150, 14138, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 28628, 0, 3,
                                                                       25988, 12248, 26138, 5150,
                                                                       5213, 14264, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 28838, 0, 3,
                                                                       26138, 12338, 26288, 5213,
                                                                       5276, 14390, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 29048, 0, 3,
                                                                       26288, 12428, 26438, 5276,
                                                                       5339, 14516, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 29258, 0, 3,
                                                                       26588, 12788, 26738, 5465,
                                                                       5528, 14894, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 29468, 0, 3,
                                                                       26738, 12878, 26888, 5528,
                                                                       5591, 15020, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 29678, 0, 3,
                                                                       26888, 12968, 27038, 5591,
                                                                       5654, 15146, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 29888, 0, 3,
                                                                       27038, 13058, 27188, 5654,
                                                                       5717, 15272, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 30098, 0, 3,
                                                                       27188, 13148, 27338, 5717,
                                                                       5780, 15398, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 30308, 0, 3,
                                                                       27338, 13238, 27488, 5780,
                                                                       5843, 15524, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 30518, 0, 3,
                                                                       27488, 13328, 27638, 5843,
                                                                       5906, 15650, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 30728, 0, 3,
                                                                       27788, 13760, 27998, 6032,
                                                                       6116, 16112, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 31008, 0, 3,
                                                                       27998, 13886, 28208, 6116,
                                                                       6200, 16280, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 31288, 0, 3,
                                                                       28208, 14012, 28418, 6200,
                                                                       6284, 16448, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 31568, 0, 3,
                                                                       28418, 14138, 28628, 6284,
                                                                       6368, 16616, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 31848, 0, 3,
                                                                       28628, 14264, 28838, 6368,
                                                                       6452, 16784, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 32128, 0, 3,
                                                                       28838, 14390, 29048, 6452,
                                                                       6536, 16952, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 32408, 0, 3,
                                                                       29258, 14894, 29468, 6704,
                                                                       6788, 17456, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 32688, 0, 3,
                                                                       29468, 15020, 29678, 6788,
                                                                       6872, 17624, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 32968, 0, 3,
                                                                       29678, 15146, 29888, 6872,
                                                                       6956, 17792, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 33248, 0, 3,
                                                                       29888, 15272, 30098, 6956,
                                                                       7040, 17960, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 33528, 0, 3,
                                                                       30098, 15398, 30308, 7040,
                                                                       7124, 18128, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 33808, 0, 3,
                                                                       30308, 15524, 30518, 7124,
                                                                       7208, 18296, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 34088, 0, 3,
                                                                       30728, 16112, 31008, 7376,
                                                                       7484, 18896, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 34448, 0, 3,
                                                                       31008, 16280, 31288, 7484,
                                                                       7592, 19112, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 34808, 0, 3,
                                                                       31288, 16448, 31568, 7592,
                                                                       7700, 19328, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 35168, 0, 3,
                                                                       31568, 16616, 31848, 7700,
                                                                       7808, 19544, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 35528, 0, 3,
                                                                       31848, 16784, 32128, 7808,
                                                                       7916, 19760, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 35888, 0, 3,
                                                                       32408, 17456, 32688, 8132,
                                                                       8240, 20408, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 36248, 0, 3,
                                                                       32688, 17624, 32968, 8240,
                                                                       8348, 20624, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 36608, 0, 3,
                                                                       32968, 17792, 33248, 8348,
                                                                       8456, 20840, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 36968, 0, 3,
                                                                       33248, 17960, 33528, 8456,
                                                                       8564, 21056, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 37328, 0, 3,
                                                                       33528, 18128, 33808, 8564,
                                                                       8672, 21272, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37688, 3, 8888,
                                                                       8894, 21488, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37703, 3, 8894,
                                                                       8900, 21498, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37718, 3, 8900,
                                                                       8906, 21508, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37733, 3, 8906,
                                                                       8912, 21518, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37748, 3, 8912,
                                                                       8918, 21528, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37763, 3, 8918,
                                                                       8924, 21538, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37778, 3, 8924,
                                                                       8930, 21548, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37793, 3, 8930,
                                                                       8936, 21558, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37808, 3, 8936,
                                                                       8942, 21568, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37823, 3, 8942,
                                                                       8948, 21578, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37838, 3, 8948,
                                                                       8954, 21588, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37853, 3, 8954,
                                                                       8960, 21598, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37868, 3, 8972,
                                                                       8978, 21608, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37883, 3, 8978,
                                                                       8984, 21618, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37898, 3, 8984,
                                                                       8990, 21628, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37913, 3, 8990,
                                                                       8996, 21638, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37928, 3, 8996,
                                                                       9002, 21648, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37943, 3, 9002,
                                                                       9008, 21658, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37958, 3, 9008,
                                                                       9014, 21668, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37973, 3, 9014,
                                                                       9020, 21678, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37988, 3, 9020,
                                                                       9026, 21688, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 38003, 3, 9026,
                                                                       9032, 21698, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 38018, 3, 9032,
                                                                       9038, 21708, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 38033, 3, 9038,
                                                                       9044, 21718, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38048, 0, 3,
                                                                       37688, 21488, 37703, 9056,
                                                                       9074, 21728, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38093, 0, 3,
                                                                       37703, 21498, 37718, 9074,
                                                                       9092, 21758, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38138, 0, 3,
                                                                       37718, 21508, 37733, 9092,
                                                                       9110, 21788, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38183, 0, 3,
                                                                       37733, 21518, 37748, 9110,
                                                                       9128, 21818, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38228, 0, 3,
                                                                       37748, 21528, 37763, 9128,
                                                                       9146, 21848, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38273, 0, 3,
                                                                       37763, 21538, 37778, 9146,
                                                                       9164, 21878, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38318, 0, 3,
                                                                       37778, 21548, 37793, 9164,
                                                                       9182, 21908, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38363, 0, 3,
                                                                       37793, 21558, 37808, 9182,
                                                                       9200, 21938, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38408, 0, 3,
                                                                       37808, 21568, 37823, 9200,
                                                                       9218, 21968, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38453, 0, 3,
                                                                       37823, 21578, 37838, 9218,
                                                                       9236, 21998, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38498, 0, 3,
                                                                       37838, 21588, 37853, 9236,
                                                                       9254, 22028, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38543, 0, 3,
                                                                       37868, 21608, 37883, 9290,
                                                                       9308, 22058, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38588, 0, 3,
                                                                       37883, 21618, 37898, 9308,
                                                                       9326, 22088, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38633, 0, 3,
                                                                       37898, 21628, 37913, 9326,
                                                                       9344, 22118, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38678, 0, 3,
                                                                       37913, 21638, 37928, 9344,
                                                                       9362, 22148, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38723, 0, 3,
                                                                       37928, 21648, 37943, 9362,
                                                                       9380, 22178, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38768, 0, 3,
                                                                       37943, 21658, 37958, 9380,
                                                                       9398, 22208, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38813, 0, 3,
                                                                       37958, 21668, 37973, 9398,
                                                                       9416, 22238, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38858, 0, 3,
                                                                       37973, 21678, 37988, 9416,
                                                                       9434, 22268, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38903, 0, 3,
                                                                       37988, 21688, 38003, 9434,
                                                                       9452, 22298, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38948, 0, 3,
                                                                       38003, 21698, 38018, 9452,
                                                                       9470, 22328, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38993, 0, 3,
                                                                       38018, 21708, 38033, 9470,
                                                                       9488, 22358, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 39038, 0, 3,
                                                                       38048, 21728, 38093, 9524,
                                                                       9560, 22388, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 39128, 0, 3,
                                                                       38093, 21758, 38138, 9560,
                                                                       9596, 22448, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 39218, 0, 3,
                                                                       38138, 21788, 38183, 9596,
                                                                       9632, 22508, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 39308, 0, 3,
                                                                       38183, 21818, 38228, 9632,
                                                                       9668, 22568, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 39398, 0, 3,
                                                                       38228, 21848, 38273, 9668,
                                                                       9704, 22628, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 39488, 0, 3,
                                                                       38273, 21878, 38318, 9704,
                                                                       9740, 22688, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 39578, 0, 3,
                                                                       38318, 21908, 38363, 9740,
                                                                       9776, 22748, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 39668, 0, 3,
                                                                       38363, 21938, 38408, 9776,
                                                                       9812, 22808, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 39758, 0, 3,
                                                                       38408, 21968, 38453, 9812,
                                                                       9848, 22868, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 39848, 0, 3,
                                                                       38453, 21998, 38498, 9848,
                                                                       9884, 22928, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 39938, 0, 3,
                                                                       38543, 22058, 38588, 9956,
                                                                       9992, 22988, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 40028, 0, 3,
                                                                       38588, 22088, 38633, 9992,
                                                                       10028, 23048, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 40118, 0, 3,
                                                                       38633, 22118, 38678,
                                                                       10028, 10064, 23108,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 40208, 0, 3,
                                                                       38678, 22148, 38723,
                                                                       10064, 10100, 23168,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 40298, 0, 3,
                                                                       38723, 22178, 38768,
                                                                       10100, 10136, 23228,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 40388, 0, 3,
                                                                       38768, 22208, 38813,
                                                                       10136, 10172, 23288,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 40478, 0, 3,
                                                                       38813, 22238, 38858,
                                                                       10172, 10208, 23348,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 40568, 0, 3,
                                                                       38858, 22268, 38903,
                                                                       10208, 10244, 23408,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 40658, 0, 3,
                                                                       38903, 22298, 38948,
                                                                       10244, 10280, 23468,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 40748, 0, 3,
                                                                       38948, 22328, 38993,
                                                                       10280, 10316, 23528,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 40838, 0, 3,
                                                                       39038, 22388, 39128,
                                                                       10388, 10448, 23588,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 40988, 0, 3,
                                                                       39128, 22448, 39218,
                                                                       10448, 10508, 23688,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 41138, 0, 3,
                                                                       39218, 22508, 39308,
                                                                       10508, 10568, 23788,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 41288, 0, 3,
                                                                       39308, 22568, 39398,
                                                                       10568, 10628, 23888,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 41438, 0, 3,
                                                                       39398, 22628, 39488,
                                                                       10628, 10688, 23988,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 41588, 0, 3,
                                                                       39488, 22688, 39578,
                                                                       10688, 10748, 24088,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 41738, 0, 3,
                                                                       39578, 22748, 39668,
                                                                       10748, 10808, 24188,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 41888, 0, 3,
                                                                       39668, 22808, 39758,
                                                                       10808, 10868, 24288,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 42038, 0, 3,
                                                                       39758, 22868, 39848,
                                                                       10868, 10928, 24388,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 42188, 0, 3,
                                                                       39938, 22988, 40028,
                                                                       11048, 11108, 24488,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 42338, 0, 3,
                                                                       40028, 23048, 40118,
                                                                       11108, 11168, 24588,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 42488, 0, 3,
                                                                       40118, 23108, 40208,
                                                                       11168, 11228, 24688,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 42638, 0, 3,
                                                                       40208, 23168, 40298,
                                                                       11228, 11288, 24788,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 42788, 0, 3,
                                                                       40298, 23228, 40388,
                                                                       11288, 11348, 24888,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 42938, 0, 3,
                                                                       40388, 23288, 40478,
                                                                       11348, 11408, 24988,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 43088, 0, 3,
                                                                       40478, 23348, 40568,
                                                                       11408, 11468, 25088,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 43238, 0, 3,
                                                                       40568, 23408, 40658,
                                                                       11468, 11528, 25188,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 43388, 0, 3,
                                                                       40658, 23468, 40748,
                                                                       11528, 11588, 25288,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 43538, 0, 3,
                                                                       40838, 23588, 40988,
                                                                       11708, 11798, 25388,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 43763, 0, 3,
                                                                       40988, 23688, 41138,
                                                                       11798, 11888, 25538,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 43988, 0, 3,
                                                                       41138, 23788, 41288,
                                                                       11888, 11978, 25688,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 44213, 0, 3,
                                                                       41288, 23888, 41438,
                                                                       11978, 12068, 25838,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 44438, 0, 3,
                                                                       41438, 23988, 41588,
                                                                       12068, 12158, 25988,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 44663, 0, 3,
                                                                       41588, 24088, 41738,
                                                                       12158, 12248, 26138,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 44888, 0, 3,
                                                                       41738, 24188, 41888,
                                                                       12248, 12338, 26288,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 45113, 0, 3,
                                                                       41888, 24288, 42038,
                                                                       12338, 12428, 26438,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 45338, 0, 3,
                                                                       42188, 24488, 42338,
                                                                       12608, 12698, 26588,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 45563, 0, 3,
                                                                       42338, 24588, 42488,
                                                                       12698, 12788, 26738,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 45788, 0, 3,
                                                                       42488, 24688, 42638,
                                                                       12788, 12878, 26888,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 46013, 0, 3,
                                                                       42638, 24788, 42788,
                                                                       12878, 12968, 27038,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 46238, 0, 3,
                                                                       42788, 24888, 42938,
                                                                       12968, 13058, 27188,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 46463, 0, 3,
                                                                       42938, 24988, 43088,
                                                                       13058, 13148, 27338,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 46688, 0, 3,
                                                                       43088, 25088, 43238,
                                                                       13148, 13238, 27488,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 46913, 0, 3,
                                                                       43238, 25188, 43388,
                                                                       13238, 13328, 27638,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 47138, 0, 3,
                                                                       43538, 25388, 43763,
                                                                       13508, 13634, 27788,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 47453, 0, 3,
                                                                       43763, 25538, 43988,
                                                                       13634, 13760, 27998,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 47768, 0, 3,
                                                                       43988, 25688, 44213,
                                                                       13760, 13886, 28208,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 48083, 0, 3,
                                                                       44213, 25838, 44438,
                                                                       13886, 14012, 28418,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 48398, 0, 3,
                                                                       44438, 25988, 44663,
                                                                       14012, 14138, 28628,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 48713, 0, 3,
                                                                       44663, 26138, 44888,
                                                                       14138, 14264, 28838,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 49028, 0, 3,
                                                                       44888, 26288, 45113,
                                                                       14264, 14390, 29048,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 49343, 0, 3,
                                                                       45338, 26588, 45563,
                                                                       14642, 14768, 29258,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 49658, 0, 3,
                                                                       45563, 26738, 45788,
                                                                       14768, 14894, 29468,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 49973, 0, 3,
                                                                       45788, 26888, 46013,
                                                                       14894, 15020, 29678,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 50288, 0, 3,
                                                                       46013, 27038, 46238,
                                                                       15020, 15146, 29888,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 50603, 0, 3,
                                                                       46238, 27188, 46463,
                                                                       15146, 15272, 30098,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 50918, 0, 3,
                                                                       46463, 27338, 46688,
                                                                       15272, 15398, 30308,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 51233, 0, 3,
                                                                       46688, 27488, 46913,
                                                                       15398, 15524, 30518,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 51548, 0, 3,
                                                                       47138, 27788, 47453,
                                                                       15776, 15944, 30728,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 51968, 0, 3,
                                                                       47453, 27998, 47768,
                                                                       15944, 16112, 31008,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 52388, 0, 3,
                                                                       47768, 28208, 48083,
                                                                       16112, 16280, 31288,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 52808, 0, 3,
                                                                       48083, 28418, 48398,
                                                                       16280, 16448, 31568,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 53228, 0, 3,
                                                                       48398, 28628, 48713,
                                                                       16448, 16616, 31848,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 53648, 0, 3,
                                                                       48713, 28838, 49028,
                                                                       16616, 16784, 32128,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 54068, 0, 3,
                                                                       49343, 29258, 49658,
                                                                       17120, 17288, 32408,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 54488, 0, 3,
                                                                       49658, 29468, 49973,
                                                                       17288, 17456, 32688,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 54908, 0, 3,
                                                                       49973, 29678, 50288,
                                                                       17456, 17624, 32968,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 55328, 0, 3,
                                                                       50288, 29888, 50603,
                                                                       17624, 17792, 33248,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 55748, 0, 3,
                                                                       50603, 30098, 50918,
                                                                       17792, 17960, 33528,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 56168, 0, 3,
                                                                       50918, 30308, 51233,
                                                                       17960, 18128, 33808,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 56588, 0, 3,
                                                                       51548, 30728, 51968,
                                                                       18464, 18680, 34088,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 57128, 0, 3,
                                                                       51968, 31008, 52388,
                                                                       18680, 18896, 34448,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 57668, 0, 3,
                                                                       52388, 31288, 52808,
                                                                       18896, 19112, 34808,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 58208, 0, 3,
                                                                       52808, 31568, 53228,
                                                                       19112, 19328, 35168,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 58748, 0, 3,
                                                                       53228, 31848, 53648,
                                                                       19328, 19544, 35528,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 59288, 0, 3,
                                                                       54068, 32408, 54488,
                                                                       19976, 20192, 35888,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 59828, 0, 3,
                                                                       54488, 32688, 54908,
                                                                       20192, 20408, 36248,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 60368, 0, 3,
                                                                       54908, 32968, 55328,
                                                                       20408, 20624, 36608,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 60908, 0, 3,
                                                                       55328, 33248, 55748,
                                                                       20624, 20840, 36968,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 61448, 0, 3,
                                                                       55748, 33528, 56168,
                                                                       20840, 21056, 37328,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 61988, 3, 21488,
                                                                       21498, 37718, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 62009, 3, 21498,
                                                                       21508, 37733, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 62030, 3, 21508,
                                                                       21518, 37748, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 62051, 3, 21518,
                                                                       21528, 37763, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 62072, 3, 21528,
                                                                       21538, 37778, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 62093, 3, 21538,
                                                                       21548, 37793, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 62114, 3, 21548,
                                                                       21558, 37808, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 62135, 3, 21558,
                                                                       21568, 37823, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 62156, 3, 21568,
                                                                       21578, 37838, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 62177, 3, 21578,
                                                                       21588, 37853, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 62198, 3, 21608,
                                                                       21618, 37898, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 62219, 3, 21618,
                                                                       21628, 37913, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 62240, 3, 21628,
                                                                       21638, 37928, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 62261, 3, 21638,
                                                                       21648, 37943, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 62282, 3, 21648,
                                                                       21658, 37958, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 62303, 3, 21658,
                                                                       21668, 37973, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 62324, 3, 21668,
                                                                       21678, 37988, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 62345, 3, 21678,
                                                                       21688, 38003, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 62366, 3, 21688,
                                                                       21698, 38018, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 62387, 3, 21698,
                                                                       21708, 38033, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 62408, 0, 3,
                                                                       61988, 37718, 62009,
                                                                       21728, 21758, 38138,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 62471, 0, 3,
                                                                       62009, 37733, 62030,
                                                                       21758, 21788, 38183,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 62534, 0, 3,
                                                                       62030, 37748, 62051,
                                                                       21788, 21818, 38228,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 62597, 0, 3,
                                                                       62051, 37763, 62072,
                                                                       21818, 21848, 38273,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 62660, 0, 3,
                                                                       62072, 37778, 62093,
                                                                       21848, 21878, 38318,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 62723, 0, 3,
                                                                       62093, 37793, 62114,
                                                                       21878, 21908, 38363,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 62786, 0, 3,
                                                                       62114, 37808, 62135,
                                                                       21908, 21938, 38408,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 62849, 0, 3,
                                                                       62135, 37823, 62156,
                                                                       21938, 21968, 38453,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 62912, 0, 3,
                                                                       62156, 37838, 62177,
                                                                       21968, 21998, 38498,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 62975, 0, 3,
                                                                       62198, 37898, 62219,
                                                                       22058, 22088, 38633,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 63038, 0, 3,
                                                                       62219, 37913, 62240,
                                                                       22088, 22118, 38678,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 63101, 0, 3,
                                                                       62240, 37928, 62261,
                                                                       22118, 22148, 38723,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 63164, 0, 3,
                                                                       62261, 37943, 62282,
                                                                       22148, 22178, 38768,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 63227, 0, 3,
                                                                       62282, 37958, 62303,
                                                                       22178, 22208, 38813,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 63290, 0, 3,
                                                                       62303, 37973, 62324,
                                                                       22208, 22238, 38858,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 63353, 0, 3,
                                                                       62324, 37988, 62345,
                                                                       22238, 22268, 38903,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 63416, 0, 3,
                                                                       62345, 38003, 62366,
                                                                       22268, 22298, 38948,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 63479, 0, 3,
                                                                       62366, 38018, 62387,
                                                                       22298, 22328, 38993,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 63542, 0, 3,
                                                                       62408, 38138, 62471,
                                                                       22388, 22448, 39218,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 63668, 0, 3,
                                                                       62471, 38183, 62534,
                                                                       22448, 22508, 39308,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 63794, 0, 3,
                                                                       62534, 38228, 62597,
                                                                       22508, 22568, 39398,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 63920, 0, 3,
                                                                       62597, 38273, 62660,
                                                                       22568, 22628, 39488,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 64046, 0, 3,
                                                                       62660, 38318, 62723,
                                                                       22628, 22688, 39578,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 64172, 0, 3,
                                                                       62723, 38363, 62786,
                                                                       22688, 22748, 39668,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 64298, 0, 3,
                                                                       62786, 38408, 62849,
                                                                       22748, 22808, 39758,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 64424, 0, 3,
                                                                       62849, 38453, 62912,
                                                                       22808, 22868, 39848,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 64550, 0, 3,
                                                                       62975, 38633, 63038,
                                                                       22988, 23048, 40118,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 64676, 0, 3,
                                                                       63038, 38678, 63101,
                                                                       23048, 23108, 40208,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 64802, 0, 3,
                                                                       63101, 38723, 63164,
                                                                       23108, 23168, 40298,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 64928, 0, 3,
                                                                       63164, 38768, 63227,
                                                                       23168, 23228, 40388,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 65054, 0, 3,
                                                                       63227, 38813, 63290,
                                                                       23228, 23288, 40478,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 65180, 0, 3,
                                                                       63290, 38858, 63353,
                                                                       23288, 23348, 40568,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 65306, 0, 3,
                                                                       63353, 38903, 63416,
                                                                       23348, 23408, 40658,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 65432, 0, 3,
                                                                       63416, 38948, 63479,
                                                                       23408, 23468, 40748,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 65558, 0, 3,
                                                                       63542, 39218, 63668,
                                                                       23588, 23688, 41138,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 65768, 0, 3,
                                                                       63668, 39308, 63794,
                                                                       23688, 23788, 41288,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 65978, 0, 3,
                                                                       63794, 39398, 63920,
                                                                       23788, 23888, 41438,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 66188, 0, 3,
                                                                       63920, 39488, 64046,
                                                                       23888, 23988, 41588,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 66398, 0, 3,
                                                                       64046, 39578, 64172,
                                                                       23988, 24088, 41738,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 66608, 0, 3,
                                                                       64172, 39668, 64298,
                                                                       24088, 24188, 41888,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 66818, 0, 3,
                                                                       64298, 39758, 64424,
                                                                       24188, 24288, 42038,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 67028, 0, 3,
                                                                       64550, 40118, 64676,
                                                                       24488, 24588, 42488,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 67238, 0, 3,
                                                                       64676, 40208, 64802,
                                                                       24588, 24688, 42638,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 67448, 0, 3,
                                                                       64802, 40298, 64928,
                                                                       24688, 24788, 42788,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 67658, 0, 3,
                                                                       64928, 40388, 65054,
                                                                       24788, 24888, 42938,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 67868, 0, 3,
                                                                       65054, 40478, 65180,
                                                                       24888, 24988, 43088,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 68078, 0, 3,
                                                                       65180, 40568, 65306,
                                                                       24988, 25088, 43238,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 68288, 0, 3,
                                                                       65306, 40658, 65432,
                                                                       25088, 25188, 43388,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 68498, 0, 3,
                                                                       65558, 41138, 65768,
                                                                       25388, 25538, 43988,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 68813, 0, 3,
                                                                       65768, 41288, 65978,
                                                                       25538, 25688, 44213,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 69128, 0, 3,
                                                                       65978, 41438, 66188,
                                                                       25688, 25838, 44438,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 69443, 0, 3,
                                                                       66188, 41588, 66398,
                                                                       25838, 25988, 44663,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 69758, 0, 3,
                                                                       66398, 41738, 66608,
                                                                       25988, 26138, 44888,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 70073, 0, 3,
                                                                       66608, 41888, 66818,
                                                                       26138, 26288, 45113,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 70388, 0, 3,
                                                                       67028, 42488, 67238,
                                                                       26588, 26738, 45788,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 70703, 0, 3,
                                                                       67238, 42638, 67448,
                                                                       26738, 26888, 46013,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 71018, 0, 3,
                                                                       67448, 42788, 67658,
                                                                       26888, 27038, 46238,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 71333, 0, 3,
                                                                       67658, 42938, 67868,
                                                                       27038, 27188, 46463,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 71648, 0, 3,
                                                                       67868, 43088, 68078,
                                                                       27188, 27338, 46688,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 71963, 0, 3,
                                                                       68078, 43238, 68288,
                                                                       27338, 27488, 46913,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 72278, 0, 3,
                                                                       68498, 43988, 68813,
                                                                       27788, 27998, 47768,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 72719, 0, 3,
                                                                       68813, 44213, 69128,
                                                                       27998, 28208, 48083,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 73160, 0, 3,
                                                                       69128, 44438, 69443,
                                                                       28208, 28418, 48398,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 73601, 0, 3,
                                                                       69443, 44663, 69758,
                                                                       28418, 28628, 48713,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 74042, 0, 3,
                                                                       69758, 44888, 70073,
                                                                       28628, 28838, 49028,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 74483, 0, 3,
                                                                       70388, 45788, 70703,
                                                                       29258, 29468, 49973,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 74924, 0, 3,
                                                                       70703, 46013, 71018,
                                                                       29468, 29678, 50288,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 75365, 0, 3,
                                                                       71018, 46238, 71333,
                                                                       29678, 29888, 50603,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 75806, 0, 3,
                                                                       71333, 46463, 71648,
                                                                       29888, 30098, 50918,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 76247, 0, 3,
                                                                       71648, 46688, 71963,
                                                                       30098, 30308, 51233,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 76688, 0, 3,
                                                                       72278, 47768, 72719,
                                                                       30728, 31008, 52388,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 77276, 0, 3,
                                                                       72719, 48083, 73160,
                                                                       31008, 31288, 52808,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 77864, 0, 3,
                                                                       73160, 48398, 73601,
                                                                       31288, 31568, 53228,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 78452, 0, 3,
                                                                       73601, 48713, 74042,
                                                                       31568, 31848, 53648,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 79040, 0, 3,
                                                                       74483, 49973, 74924,
                                                                       32408, 32688, 54908,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 79628, 0, 3,
                                                                       74924, 50288, 75365,
                                                                       32688, 32968, 55328,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 80216, 0, 3,
                                                                       75365, 50603, 75806,
                                                                       32968, 33248, 55748,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 80804, 0, 3,
                                                                       75806, 50918, 76247,
                                                                       33248, 33528, 56168,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 81392, 0, 3,
                                                                       76688, 52388, 77276,
                                                                       34088, 34448, 57668,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 82148, 0, 3,
                                                                       77276, 52808, 77864,
                                                                       34448, 34808, 58208,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 82904, 0, 3,
                                                                       77864, 53228, 78452,
                                                                       34808, 35168, 58748,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 83660, 0, 3,
                                                                       79040, 54908, 79628,
                                                                       35888, 36248, 60368,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 84416, 0, 3,
                                                                       79628, 55328, 80216,
                                                                       36248, 36608, 60908,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 85172, 0, 3,
                                                                       80216, 55748, 80804,
                                                                       36608, 36968, 61448,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 85928, 3, 37688,
                                                                       37703, 61988, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 85956, 3, 37703,
                                                                       37718, 62009, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 85984, 3, 37718,
                                                                       37733, 62030, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 86012, 3, 37733,
                                                                       37748, 62051, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 86040, 3, 37748,
                                                                       37763, 62072, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 86068, 3, 37763,
                                                                       37778, 62093, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 86096, 3, 37778,
                                                                       37793, 62114, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 86124, 3, 37793,
                                                                       37808, 62135, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 86152, 3, 37808,
                                                                       37823, 62156, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 86180, 3, 37823,
                                                                       37838, 62177, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 86208, 3, 37868,
                                                                       37883, 62198, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 86236, 3, 37883,
                                                                       37898, 62219, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 86264, 3, 37898,
                                                                       37913, 62240, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 86292, 3, 37913,
                                                                       37928, 62261, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 86320, 3, 37928,
                                                                       37943, 62282, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 86348, 3, 37943,
                                                                       37958, 62303, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 86376, 3, 37958,
                                                                       37973, 62324, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 86404, 3, 37973,
                                                                       37988, 62345, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 86432, 3, 37988,
                                                                       38003, 62366, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 86460, 3, 38003,
                                                                       38018, 62387, ncols,
                                                                       gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 86488, 0, 3,
                                                                       85928, 61988, 85956,
                                                                       38048, 38093, 62408,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 86572, 0, 3,
                                                                       85956, 62009, 85984,
                                                                       38093, 38138, 62471,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 86656, 0, 3,
                                                                       85984, 62030, 86012,
                                                                       38138, 38183, 62534,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 86740, 0, 3,
                                                                       86012, 62051, 86040,
                                                                       38183, 38228, 62597,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 86824, 0, 3,
                                                                       86040, 62072, 86068,
                                                                       38228, 38273, 62660,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 86908, 0, 3,
                                                                       86068, 62093, 86096,
                                                                       38273, 38318, 62723,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 86992, 0, 3,
                                                                       86096, 62114, 86124,
                                                                       38318, 38363, 62786,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 87076, 0, 3,
                                                                       86124, 62135, 86152,
                                                                       38363, 38408, 62849,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 87160, 0, 3,
                                                                       86152, 62156, 86180,
                                                                       38408, 38453, 62912,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 87244, 0, 3,
                                                                       86208, 62198, 86236,
                                                                       38543, 38588, 62975,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 87328, 0, 3,
                                                                       86236, 62219, 86264,
                                                                       38588, 38633, 63038,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 87412, 0, 3,
                                                                       86264, 62240, 86292,
                                                                       38633, 38678, 63101,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 87496, 0, 3,
                                                                       86292, 62261, 86320,
                                                                       38678, 38723, 63164,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 87580, 0, 3,
                                                                       86320, 62282, 86348,
                                                                       38723, 38768, 63227,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 87664, 0, 3,
                                                                       86348, 62303, 86376,
                                                                       38768, 38813, 63290,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 87748, 0, 3,
                                                                       86376, 62324, 86404,
                                                                       38813, 38858, 63353,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 87832, 0, 3,
                                                                       86404, 62345, 86432,
                                                                       38858, 38903, 63416,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 87916, 0, 3,
                                                                       86432, 62366, 86460,
                                                                       38903, 38948, 63479,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 88000, 0, 3,
                                                                       86488, 62408, 86572,
                                                                       39038, 39128, 63542,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 88168, 0, 3,
                                                                       86572, 62471, 86656,
                                                                       39128, 39218, 63668,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 88336, 0, 3,
                                                                       86656, 62534, 86740,
                                                                       39218, 39308, 63794,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 88504, 0, 3,
                                                                       86740, 62597, 86824,
                                                                       39308, 39398, 63920,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 88672, 0, 3,
                                                                       86824, 62660, 86908,
                                                                       39398, 39488, 64046,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 88840, 0, 3,
                                                                       86908, 62723, 86992,
                                                                       39488, 39578, 64172,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 89008, 0, 3,
                                                                       86992, 62786, 87076,
                                                                       39578, 39668, 64298,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 89176, 0, 3,
                                                                       87076, 62849, 87160,
                                                                       39668, 39758, 64424,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 89344, 0, 3,
                                                                       87244, 62975, 87328,
                                                                       39938, 40028, 64550,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 89512, 0, 3,
                                                                       87328, 63038, 87412,
                                                                       40028, 40118, 64676,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 89680, 0, 3,
                                                                       87412, 63101, 87496,
                                                                       40118, 40208, 64802,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 89848, 0, 3,
                                                                       87496, 63164, 87580,
                                                                       40208, 40298, 64928,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 90016, 0, 3,
                                                                       87580, 63227, 87664,
                                                                       40298, 40388, 65054,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 90184, 0, 3,
                                                                       87664, 63290, 87748,
                                                                       40388, 40478, 65180,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 90352, 0, 3,
                                                                       87748, 63353, 87832,
                                                                       40478, 40568, 65306,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 90520, 0, 3,
                                                                       87832, 63416, 87916,
                                                                       40568, 40658, 65432,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 90688, 0, 3,
                                                                       88000, 63542, 88168,
                                                                       40838, 40988, 65558,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 90968, 0, 3,
                                                                       88168, 63668, 88336,
                                                                       40988, 41138, 65768,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 91248, 0, 3,
                                                                       88336, 63794, 88504,
                                                                       41138, 41288, 65978,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 91528, 0, 3,
                                                                       88504, 63920, 88672,
                                                                       41288, 41438, 66188,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 91808, 0, 3,
                                                                       88672, 64046, 88840,
                                                                       41438, 41588, 66398,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 92088, 0, 3,
                                                                       88840, 64172, 89008,
                                                                       41588, 41738, 66608,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 92368, 0, 3,
                                                                       89008, 64298, 89176,
                                                                       41738, 41888, 66818,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 92648, 0, 3,
                                                                       89344, 64550, 89512,
                                                                       42188, 42338, 67028,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 92928, 0, 3,
                                                                       89512, 64676, 89680,
                                                                       42338, 42488, 67238,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 93208, 0, 3,
                                                                       89680, 64802, 89848,
                                                                       42488, 42638, 67448,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 93488, 0, 3,
                                                                       89848, 64928, 90016,
                                                                       42638, 42788, 67658,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 93768, 0, 3,
                                                                       90016, 65054, 90184,
                                                                       42788, 42938, 67868,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 94048, 0, 3,
                                                                       90184, 65180, 90352,
                                                                       42938, 43088, 68078,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 94328, 0, 3,
                                                                       90352, 65306, 90520,
                                                                       43088, 43238, 68288,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 94608, 0, 3,
                                                                       90688, 65558, 90968,
                                                                       43538, 43763, 68498,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 95028, 0, 3,
                                                                       90968, 65768, 91248,
                                                                       43763, 43988, 68813,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 95448, 0, 3,
                                                                       91248, 65978, 91528,
                                                                       43988, 44213, 69128,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 95868, 0, 3,
                                                                       91528, 66188, 91808,
                                                                       44213, 44438, 69443,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 96288, 0, 3,
                                                                       91808, 66398, 92088,
                                                                       44438, 44663, 69758,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 96708, 0, 3,
                                                                       92088, 66608, 92368,
                                                                       44663, 44888, 70073,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 97128, 0, 3,
                                                                       92648, 67028, 92928,
                                                                       45338, 45563, 70388,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 97548, 0, 3,
                                                                       92928, 67238, 93208,
                                                                       45563, 45788, 70703,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 97968, 0, 3,
                                                                       93208, 67448, 93488,
                                                                       45788, 46013, 71018,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 98388, 0, 3,
                                                                       93488, 67658, 93768,
                                                                       46013, 46238, 71333,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 98808, 0, 3,
                                                                       93768, 67868, 94048,
                                                                       46238, 46463, 71648,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 99228, 0, 3,
                                                                       94048, 68078, 94328,
                                                                       46463, 46688, 71963,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 99648, 0, 3,
                                                                       94608, 68498, 95028,
                                                                       47138, 47453, 72278,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 100236, 0, 3,
                                                                       95028, 68813, 95448,
                                                                       47453, 47768, 72719,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 100824, 0, 3,
                                                                       95448, 69128, 95868,
                                                                       47768, 48083, 73160,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 101412, 0, 3,
                                                                       95868, 69443, 96288,
                                                                       48083, 48398, 73601,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 102000, 0, 3,
                                                                       96288, 69758, 96708,
                                                                       48398, 48713, 74042,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 102588, 0, 3,
                                                                       97128, 70388, 97548,
                                                                       49343, 49658, 74483,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 103176, 0, 3,
                                                                       97548, 70703, 97968,
                                                                       49658, 49973, 74924,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 103764, 0, 3,
                                                                       97968, 71018, 98388,
                                                                       49973, 50288, 75365,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 104352, 0, 3,
                                                                       98388, 71333, 98808,
                                                                       50288, 50603, 75806,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 104940, 0, 3,
                                                                       98808, 71648, 99228,
                                                                       50603, 50918, 76247,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 105528, 0, 3,
                                                                       99648, 72278, 100236,
                                                                       51548, 51968, 76688,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 106312, 0, 3,
                                                                       100236, 72719, 100824,
                                                                       51968, 52388, 77276,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 107096, 0, 3,
                                                                       100824, 73160, 101412,
                                                                       52388, 52808, 77864,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 107880, 0, 3,
                                                                       101412, 73601, 102000,
                                                                       52808, 53228, 78452,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 108664, 0, 3,
                                                                       102588, 74483, 103176,
                                                                       54068, 54488, 79040,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 109448, 0, 3,
                                                                       103176, 74924, 103764,
                                                                       54488, 54908, 79628,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 110232, 0, 3,
                                                                       103764, 75365, 104352,
                                                                       54908, 55328, 80216,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 111016, 0, 3,
                                                                       104352, 75806, 104940,
                                                                       55328, 55748, 80804,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 111800, 0, 3,
                                                                       105528, 76688, 106312,
                                                                       56588, 57128, 81392,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 112808, 0, 3,
                                                                       106312, 77276, 107096,
                                                                       57128, 57668, 82148,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 113816, 0, 3,
                                                                       107096, 77864, 107880,
                                                                       57668, 58208, 82904,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 114824, 0, 3,
                                                                       108664, 79040, 109448,
                                                                       59288, 59828, 83660,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 115832, 0, 3,
                                                                       109448, 79628, 110232,
                                                                       59828, 60368, 84416,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 116840, 0, 3,
                                                                       110232, 80216, 111016,
                                                                       60368, 60908, 85172,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 117848, 3, 61988,
                                                                       62009, 85984, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 117884, 3, 62009,
                                                                       62030, 86012, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 117920, 3, 62030,
                                                                       62051, 86040, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 117956, 3, 62051,
                                                                       62072, 86068, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 117992, 3, 62072,
                                                                       62093, 86096, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 118028, 3, 62093,
                                                                       62114, 86124, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 118064, 3, 62114,
                                                                       62135, 86152, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 118100, 3, 62135,
                                                                       62156, 86180, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 118136, 3, 62198,
                                                                       62219, 86264, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 118172, 3, 62219,
                                                                       62240, 86292, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 118208, 3, 62240,
                                                                       62261, 86320, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 118244, 3, 62261,
                                                                       62282, 86348, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 118280, 3, 62282,
                                                                       62303, 86376, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 118316, 3, 62303,
                                                                       62324, 86404, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 118352, 3, 62324,
                                                                       62345, 86432, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 118388, 3, 62345,
                                                                       62366, 86460, ncols,
                                                                       gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 118424, 0, 3,
                                                                       117848, 85984, 117884,
                                                                       62408, 62471, 86656,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 118532, 0, 3,
                                                                       117884, 86012, 117920,
                                                                       62471, 62534, 86740,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 118640, 0, 3,
                                                                       117920, 86040, 117956,
                                                                       62534, 62597, 86824,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 118748, 0, 3,
                                                                       117956, 86068, 117992,
                                                                       62597, 62660, 86908,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 118856, 0, 3,
                                                                       117992, 86096, 118028,
                                                                       62660, 62723, 86992,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 118964, 0, 3,
                                                                       118028, 86124, 118064,
                                                                       62723, 62786, 87076,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 119072, 0, 3,
                                                                       118064, 86152, 118100,
                                                                       62786, 62849, 87160,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 119180, 0, 3,
                                                                       118136, 86264, 118172,
                                                                       62975, 63038, 87412,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 119288, 0, 3,
                                                                       118172, 86292, 118208,
                                                                       63038, 63101, 87496,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 119396, 0, 3,
                                                                       118208, 86320, 118244,
                                                                       63101, 63164, 87580,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 119504, 0, 3,
                                                                       118244, 86348, 118280,
                                                                       63164, 63227, 87664,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 119612, 0, 3,
                                                                       118280, 86376, 118316,
                                                                       63227, 63290, 87748,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 119720, 0, 3,
                                                                       118316, 86404, 118352,
                                                                       63290, 63353, 87832,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 119828, 0, 3,
                                                                       118352, 86432, 118388,
                                                                       63353, 63416, 87916,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 119936, 0, 3,
                                                                       118424, 86656, 118532,
                                                                       63542, 63668, 88336,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 120152, 0, 3,
                                                                       118532, 86740, 118640,
                                                                       63668, 63794, 88504,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 120368, 0, 3,
                                                                       118640, 86824, 118748,
                                                                       63794, 63920, 88672,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 120584, 0, 3,
                                                                       118748, 86908, 118856,
                                                                       63920, 64046, 88840,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 120800, 0, 3,
                                                                       118856, 86992, 118964,
                                                                       64046, 64172, 89008,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 121016, 0, 3,
                                                                       118964, 87076, 119072,
                                                                       64172, 64298, 89176,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 121232, 0, 3,
                                                                       119180, 87412, 119288,
                                                                       64550, 64676, 89680,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 121448, 0, 3,
                                                                       119288, 87496, 119396,
                                                                       64676, 64802, 89848,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 121664, 0, 3,
                                                                       119396, 87580, 119504,
                                                                       64802, 64928, 90016,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 121880, 0, 3,
                                                                       119504, 87664, 119612,
                                                                       64928, 65054, 90184,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 122096, 0, 3,
                                                                       119612, 87748, 119720,
                                                                       65054, 65180, 90352,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 122312, 0, 3,
                                                                       119720, 87832, 119828,
                                                                       65180, 65306, 90520,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 122528, 0, 3,
                                                                       119936, 88336, 120152,
                                                                       65558, 65768, 91248,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 122888, 0, 3,
                                                                       120152, 88504, 120368,
                                                                       65768, 65978, 91528,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 123248, 0, 3,
                                                                       120368, 88672, 120584,
                                                                       65978, 66188, 91808,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 123608, 0, 3,
                                                                       120584, 88840, 120800,
                                                                       66188, 66398, 92088,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 123968, 0, 3,
                                                                       120800, 89008, 121016,
                                                                       66398, 66608, 92368,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 124328, 0, 3,
                                                                       121232, 89680, 121448,
                                                                       67028, 67238, 93208,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 124688, 0, 3,
                                                                       121448, 89848, 121664,
                                                                       67238, 67448, 93488,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 125048, 0, 3,
                                                                       121664, 90016, 121880,
                                                                       67448, 67658, 93768,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 125408, 0, 3,
                                                                       121880, 90184, 122096,
                                                                       67658, 67868, 94048,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 125768, 0, 3,
                                                                       122096, 90352, 122312,
                                                                       67868, 68078, 94328,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 126128, 0, 3,
                                                                       122528, 91248, 122888,
                                                                       68498, 68813, 95448,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 126668, 0, 3,
                                                                       122888, 91528, 123248,
                                                                       68813, 69128, 95868,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 127208, 0, 3,
                                                                       123248, 91808, 123608,
                                                                       69128, 69443, 96288,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 127748, 0, 3,
                                                                       123608, 92088, 123968,
                                                                       69443, 69758, 96708,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 128288, 0, 3,
                                                                       124328, 93208, 124688,
                                                                       70388, 70703, 97968,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 128828, 0, 3,
                                                                       124688, 93488, 125048,
                                                                       70703, 71018, 98388,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 129368, 0, 3,
                                                                       125048, 93768, 125408,
                                                                       71018, 71333, 98808,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 129908, 0, 3,
                                                                       125408, 94048, 125768,
                                                                       71333, 71648, 99228,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 130448, 0, 3,
                                                                       126128, 95448, 126668,
                                                                       72278, 72719, 100824,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 131204, 0, 3,
                                                                       126668, 95868, 127208,
                                                                       72719, 73160, 101412,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 131960, 0, 3,
                                                                       127208, 96288, 127748,
                                                                       73160, 73601, 102000,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 132716, 0, 3,
                                                                       128288, 97968, 128828,
                                                                       74483, 74924, 103764,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 133472, 0, 3,
                                                                       128828, 98388, 129368,
                                                                       74924, 75365, 104352,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 134228, 0, 3,
                                                                       129368, 98808, 129908,
                                                                       75365, 75806, 104940,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 134984, 0, 3,
                                                                       130448, 100824, 131204,
                                                                       76688, 77276, 107096,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 135992, 0, 3,
                                                                       131204, 101412, 131960,
                                                                       77276, 77864, 107880,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 137000, 0, 3,
                                                                       132716, 103764, 133472,
                                                                       79040, 79628, 110232,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 138008, 0, 3,
                                                                       133472, 104352, 134228,
                                                                       79628, 80216, 111016,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 139016, 0, 3,
                                                                       134984, 107096, 135992,
                                                                       81392, 82148, 113816,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 140312, 0, 3,
                                                                       137000, 110232, 138008,
                                                                       83660, 84416, 116840,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 141608, 3, 85928,
                                                                       85956, 117848, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 141653, 3, 85956,
                                                                       85984, 117884, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 141698, 3, 85984,
                                                                       86012, 117920, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 141743, 3, 86012,
                                                                       86040, 117956, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 141788, 3, 86040,
                                                                       86068, 117992, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 141833, 3, 86068,
                                                                       86096, 118028, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 141878, 3, 86096,
                                                                       86124, 118064, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 141923, 3, 86124,
                                                                       86152, 118100, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 141968, 3, 86208,
                                                                       86236, 118136, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 142013, 3, 86236,
                                                                       86264, 118172, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 142058, 3, 86264,
                                                                       86292, 118208, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 142103, 3, 86292,
                                                                       86320, 118244, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 142148, 3, 86320,
                                                                       86348, 118280, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 142193, 3, 86348,
                                                                       86376, 118316, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 142238, 3, 86376,
                                                                       86404, 118352, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 142283, 3, 86404,
                                                                       86432, 118388, ncols,
                                                                       gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 142328, 0, 3,
                                                                       141608, 117848, 141653,
                                                                       86488, 86572, 118424,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 142463, 0, 3,
                                                                       141653, 117884, 141698,
                                                                       86572, 86656, 118532,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 142598, 0, 3,
                                                                       141698, 117920, 141743,
                                                                       86656, 86740, 118640,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 142733, 0, 3,
                                                                       141743, 117956, 141788,
                                                                       86740, 86824, 118748,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 142868, 0, 3,
                                                                       141788, 117992, 141833,
                                                                       86824, 86908, 118856,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 143003, 0, 3,
                                                                       141833, 118028, 141878,
                                                                       86908, 86992, 118964,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 143138, 0, 3,
                                                                       141878, 118064, 141923,
                                                                       86992, 87076, 119072,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 143273, 0, 3,
                                                                       141968, 118136, 142013,
                                                                       87244, 87328, 119180,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 143408, 0, 3,
                                                                       142013, 118172, 142058,
                                                                       87328, 87412, 119288,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 143543, 0, 3,
                                                                       142058, 118208, 142103,
                                                                       87412, 87496, 119396,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 143678, 0, 3,
                                                                       142103, 118244, 142148,
                                                                       87496, 87580, 119504,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 143813, 0, 3,
                                                                       142148, 118280, 142193,
                                                                       87580, 87664, 119612,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 143948, 0, 3,
                                                                       142193, 118316, 142238,
                                                                       87664, 87748, 119720,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 144083, 0, 3,
                                                                       142238, 118352, 142283,
                                                                       87748, 87832, 119828,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 144218, 0, 3,
                                                                       142328, 118424, 142463,
                                                                       88000, 88168, 119936,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 144488, 0, 3,
                                                                       142463, 118532, 142598,
                                                                       88168, 88336, 120152,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 144758, 0, 3,
                                                                       142598, 118640, 142733,
                                                                       88336, 88504, 120368,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 145028, 0, 3,
                                                                       142733, 118748, 142868,
                                                                       88504, 88672, 120584,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 145298, 0, 3,
                                                                       142868, 118856, 143003,
                                                                       88672, 88840, 120800,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 145568, 0, 3,
                                                                       143003, 118964, 143138,
                                                                       88840, 89008, 121016,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 145838, 0, 3,
                                                                       143273, 119180, 143408,
                                                                       89344, 89512, 121232,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 146108, 0, 3,
                                                                       143408, 119288, 143543,
                                                                       89512, 89680, 121448,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 146378, 0, 3,
                                                                       143543, 119396, 143678,
                                                                       89680, 89848, 121664,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 146648, 0, 3,
                                                                       143678, 119504, 143813,
                                                                       89848, 90016, 121880,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 146918, 0, 3,
                                                                       143813, 119612, 143948,
                                                                       90016, 90184, 122096,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 147188, 0, 3,
                                                                       143948, 119720, 144083,
                                                                       90184, 90352, 122312,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 147458, 0, 3,
                                                                       144218, 119936, 144488,
                                                                       90688, 90968, 122528,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 147908, 0, 3,
                                                                       144488, 120152, 144758,
                                                                       90968, 91248, 122888,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 148358, 0, 3,
                                                                       144758, 120368, 145028,
                                                                       91248, 91528, 123248,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 148808, 0, 3,
                                                                       145028, 120584, 145298,
                                                                       91528, 91808, 123608,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 149258, 0, 3,
                                                                       145298, 120800, 145568,
                                                                       91808, 92088, 123968,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 149708, 0, 3,
                                                                       145838, 121232, 146108,
                                                                       92648, 92928, 124328,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 150158, 0, 3,
                                                                       146108, 121448, 146378,
                                                                       92928, 93208, 124688,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 150608, 0, 3,
                                                                       146378, 121664, 146648,
                                                                       93208, 93488, 125048,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 151058, 0, 3,
                                                                       146648, 121880, 146918,
                                                                       93488, 93768, 125408,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 151508, 0, 3,
                                                                       146918, 122096, 147188,
                                                                       93768, 94048, 125768,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 151958, 0, 3,
                                                                       147458, 122528, 147908,
                                                                       94608, 95028, 126128,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 152633, 0, 3,
                                                                       147908, 122888, 148358,
                                                                       95028, 95448, 126668,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 153308, 0, 3,
                                                                       148358, 123248, 148808,
                                                                       95448, 95868, 127208,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 153983, 0, 3,
                                                                       148808, 123608, 149258,
                                                                       95868, 96288, 127748,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 154658, 0, 3,
                                                                       149708, 124328, 150158,
                                                                       97128, 97548, 128288,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 155333, 0, 3,
                                                                       150158, 124688, 150608,
                                                                       97548, 97968, 128828,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 156008, 0, 3,
                                                                       150608, 125048, 151058,
                                                                       97968, 98388, 129368,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 156683, 0, 3,
                                                                       151058, 125408, 151508,
                                                                       98388, 98808, 129908,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 157358, 0, 3,
                                                                       151958, 126128, 152633,
                                                                       99648, 100236, 130448,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 158303, 0, 3,
                                                                       152633, 126668, 153308,
                                                                       100236, 100824, 131204,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 159248, 0, 3,
                                                                       153308, 127208, 153983,
                                                                       100824, 101412, 131960,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 160193, 0, 3,
                                                                       154658, 128288, 155333,
                                                                       102588, 103176, 132716,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 161138, 0, 3,
                                                                       155333, 128828, 156008,
                                                                       103176, 103764, 133472,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 162083, 0, 3,
                                                                       156008, 129368, 156683,
                                                                       103764, 104352, 134228,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 163028, 0, 3,
                                                                       157358, 130448, 158303,
                                                                       105528, 106312, 134984,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 164288, 0, 3,
                                                                       158303, 131204, 159248,
                                                                       106312, 107096, 135992,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 165548, 0, 3,
                                                                       160193, 132716, 161138,
                                                                       108664, 109448, 137000,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 166808, 0, 3,
                                                                       161138, 133472, 162083,
                                                                       109448, 110232, 138008,
                                                                       ncols, gamma, p, q);

                    compute_prim_skl_three_center_electron_repulsion_0(buffer, 168068, 0, 3,
                                                                       163028, 134984, 164288,
                                                                       111800, 112808, 139016,
                                                                       ncols, gamma, p, q);

                    compute_prim_skl_three_center_electron_repulsion_0(buffer, 169688, 0, 3,
                                                                       165548, 137000, 166808,
                                                                       114824, 115832, 140312,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 171308, 151958, 675, ncols);

                    simdfunc::contract_primitives(buffer, 172238, 154658, 675, ncols);

                    simdfunc::contract_primitives(buffer, 173168, 157358, 945, ncols);

                    simdfunc::contract_primitives(buffer, 174470, 160193, 945, ncols);

                    simdfunc::contract_primitives(buffer, 175772, 163028, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 177508, 165548, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 179244, 168068, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 181476, 169688, 1620, ncols);
                }
            }
        }

        simdtrf::transform_l_inner(buffer, 171983, 171308, 15, 1, nmax);

        simdtrf::transform_l_inner(buffer, 172913, 172238, 15, 1, nmax);

        simdtrf::transform_l_inner(buffer, 174113, 173168, 21, 1, nmax);

        simdtrf::transform_l_inner(buffer, 175415, 174470, 21, 1, nmax);

        simdtrf::transform_l_inner(buffer, 177032, 175772, 28, 1, nmax);

        simdtrf::transform_l_inner(buffer, 178768, 177508, 28, 1, nmax);

        simdtrf::transform_l_inner(buffer, 180864, 179244, 36, 1, nmax);

        simdtrf::transform_l_inner(buffer, 183096, 181476, 36, 1, nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 183708, 171983, 174113, 17, nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 184473, 172913, 175415, 17, nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 185238, 174113, 177032, 17, nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 186309, 175415, 178768, 17, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 187380, 177032, 180864, 17, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 188808, 178768, 183096, 17, nmax);

        simdtrf::compute_hrr_dg(buffer, coordinates, 190236, 183708, 185238, 17, nmax);

        simdtrf::compute_hrr_dg(buffer, coordinates, 191766, 184473, 186309, 17, nmax);

        simdtrf::compute_hrr_dh(buffer, coordinates, 193296, 185238, 187380, 17, nmax);

        simdtrf::compute_hrr_dh(buffer, coordinates, 195438, 186309, 188808, 17, nmax);

        simdtrf::compute_hrr_fg(buffer, coordinates, 197580, 190236, 193296, 17, nmax);

        simdtrf::compute_hrr_fg(buffer, coordinates, 200130, 191766, 195438, 17, nmax);

        simdtrf::transform_g_inner(buffer, 202680, 200130, 10, 17, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 202680, 153, nmax);

        simdtrf::transform_g_inner(buffer, 202680, 197580, 10, 17, nmax);

        simdtrf::transform_f_outer(values + 1071 * nvalues + n * npairs, nvalues, buffer, 202680,
                                   153, nmax);
    }

    for (size_t m = 0; m < 2142; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
