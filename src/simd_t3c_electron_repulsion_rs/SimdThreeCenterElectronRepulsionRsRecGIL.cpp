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


#include "SimdThreeCenterElectronRepulsionRsRecGIL.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSLD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSML.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSND.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNS.hpp"
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
#include "SimdTransformI.hpp"
#include "SimdTransformL.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_gil_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_gil_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 585897, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 3978 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 585897, 479058, 27398, dimensions);

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

                    simdfunc::compute_full_t3c_erf_boys_function(buffer, coordinates, 6, 3, 18,
                                                                 ncols, fj, i * nprim_b + j, fq,
                                                                 omega);

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 26, 3, 18,
                                                             ncols, fj, i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 46, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 49, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 52, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 55, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 58, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 61, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 64, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 67, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 70, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 73, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 76, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 79, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 82, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 85, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 88, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 91, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 94, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 97, 0, 3, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 100, 0, 3, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 103, 0, 3, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 106, 0, 3, 29, 30,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 109, 0, 3, 30, 31,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 112, 0, 3, 31, 32,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 115, 0, 3, 32, 33,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 118, 0, 3, 33, 34,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 121, 0, 3, 34, 35,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 124, 0, 3, 35, 36,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 127, 0, 3, 36, 37,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 130, 0, 3, 37, 38,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 133, 0, 3, 38, 39,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 136, 0, 3, 39, 40,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 139, 0, 3, 40, 41,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 142, 0, 3, 41, 42,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 145, 0, 3, 42, 43,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 148, 0, 3, 43, 44,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 151, 0, 3, 44, 45,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 154, 0, 3, 7, 8,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 160, 0, 3, 8, 9,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 166, 0, 3, 9, 10,
                                                                       52, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 172, 0, 3, 10, 11,
                                                                       55, 58, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 178, 0, 3, 11, 12,
                                                                       58, 61, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 184, 0, 3, 12, 13,
                                                                       61, 64, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 190, 0, 3, 13, 14,
                                                                       64, 67, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 196, 0, 3, 14, 15,
                                                                       67, 70, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 202, 0, 3, 15, 16,
                                                                       70, 73, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 208, 0, 3, 16, 17,
                                                                       73, 76, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 214, 0, 3, 17, 18,
                                                                       76, 79, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 220, 0, 3, 18, 19,
                                                                       79, 82, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 226, 0, 3, 19, 20,
                                                                       82, 85, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 232, 0, 3, 20, 21,
                                                                       85, 88, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 238, 0, 3, 21, 22,
                                                                       88, 91, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 244, 0, 3, 22, 23,
                                                                       91, 94, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 250, 0, 3, 23, 24,
                                                                       94, 97, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 256, 0, 3, 27, 28,
                                                                       100, 103, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 262, 0, 3, 28, 29,
                                                                       103, 106, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 268, 0, 3, 29, 30,
                                                                       106, 109, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 274, 0, 3, 30, 31,
                                                                       109, 112, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 280, 0, 3, 31, 32,
                                                                       112, 115, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 286, 0, 3, 32, 33,
                                                                       115, 118, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 292, 0, 3, 33, 34,
                                                                       118, 121, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 298, 0, 3, 34, 35,
                                                                       121, 124, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 304, 0, 3, 35, 36,
                                                                       124, 127, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 310, 0, 3, 36, 37,
                                                                       127, 130, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 316, 0, 3, 37, 38,
                                                                       130, 133, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 322, 0, 3, 38, 39,
                                                                       133, 136, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 328, 0, 3, 39, 40,
                                                                       136, 139, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 334, 0, 3, 40, 41,
                                                                       139, 142, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 340, 0, 3, 41, 42,
                                                                       142, 145, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 346, 0, 3, 42, 43,
                                                                       145, 148, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 352, 0, 3, 43, 44,
                                                                       148, 151, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 358, 0, 3, 46, 49,
                                                                       154, 160, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 368, 0, 3, 49, 52,
                                                                       160, 166, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 378, 0, 3, 52, 55,
                                                                       166, 172, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 388, 0, 3, 55, 58,
                                                                       172, 178, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 398, 0, 3, 58, 61,
                                                                       178, 184, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 408, 0, 3, 61, 64,
                                                                       184, 190, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 418, 0, 3, 64, 67,
                                                                       190, 196, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 428, 0, 3, 67, 70,
                                                                       196, 202, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 438, 0, 3, 70, 73,
                                                                       202, 208, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 448, 0, 3, 73, 76,
                                                                       208, 214, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 458, 0, 3, 76, 79,
                                                                       214, 220, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 468, 0, 3, 79, 82,
                                                                       220, 226, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 478, 0, 3, 82, 85,
                                                                       226, 232, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 488, 0, 3, 85, 88,
                                                                       232, 238, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 498, 0, 3, 88, 91,
                                                                       238, 244, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 508, 0, 3, 91, 94,
                                                                       244, 250, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 518, 0, 3, 100,
                                                                       103, 256, 262, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 528, 0, 3, 103,
                                                                       106, 262, 268, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 538, 0, 3, 106,
                                                                       109, 268, 274, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 548, 0, 3, 109,
                                                                       112, 274, 280, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 558, 0, 3, 112,
                                                                       115, 280, 286, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 568, 0, 3, 115,
                                                                       118, 286, 292, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 578, 0, 3, 118,
                                                                       121, 292, 298, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 588, 0, 3, 121,
                                                                       124, 298, 304, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 598, 0, 3, 124,
                                                                       127, 304, 310, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 608, 0, 3, 127,
                                                                       130, 310, 316, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 618, 0, 3, 130,
                                                                       133, 316, 322, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 628, 0, 3, 133,
                                                                       136, 322, 328, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 638, 0, 3, 136,
                                                                       139, 328, 334, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 648, 0, 3, 139,
                                                                       142, 334, 340, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 658, 0, 3, 142,
                                                                       145, 340, 346, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 668, 0, 3, 145,
                                                                       148, 346, 352, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 678, 0, 3, 154,
                                                                       160, 358, 368, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 693, 0, 3, 160,
                                                                       166, 368, 378, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 708, 0, 3, 166,
                                                                       172, 378, 388, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 723, 0, 3, 172,
                                                                       178, 388, 398, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 738, 0, 3, 178,
                                                                       184, 398, 408, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 753, 0, 3, 184,
                                                                       190, 408, 418, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 768, 0, 3, 190,
                                                                       196, 418, 428, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 783, 0, 3, 196,
                                                                       202, 428, 438, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 798, 0, 3, 202,
                                                                       208, 438, 448, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 813, 0, 3, 208,
                                                                       214, 448, 458, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 828, 0, 3, 214,
                                                                       220, 458, 468, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 843, 0, 3, 220,
                                                                       226, 468, 478, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 858, 0, 3, 226,
                                                                       232, 478, 488, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 873, 0, 3, 232,
                                                                       238, 488, 498, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 888, 0, 3, 238,
                                                                       244, 498, 508, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 903, 0, 3, 256,
                                                                       262, 518, 528, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 918, 0, 3, 262,
                                                                       268, 528, 538, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 933, 0, 3, 268,
                                                                       274, 538, 548, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 948, 0, 3, 274,
                                                                       280, 548, 558, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 963, 0, 3, 280,
                                                                       286, 558, 568, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 978, 0, 3, 286,
                                                                       292, 568, 578, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 993, 0, 3, 292,
                                                                       298, 578, 588, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 1008, 0, 3, 298,
                                                                       304, 588, 598, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 1023, 0, 3, 304,
                                                                       310, 598, 608, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 1038, 0, 3, 310,
                                                                       316, 608, 618, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 1053, 0, 3, 316,
                                                                       322, 618, 628, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 1068, 0, 3, 322,
                                                                       328, 628, 638, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 1083, 0, 3, 328,
                                                                       334, 638, 648, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 1098, 0, 3, 334,
                                                                       340, 648, 658, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 1113, 0, 3, 340,
                                                                       346, 658, 668, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1128, 0, 3, 358,
                                                                       368, 678, 693, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1149, 0, 3, 368,
                                                                       378, 693, 708, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1170, 0, 3, 378,
                                                                       388, 708, 723, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1191, 0, 3, 388,
                                                                       398, 723, 738, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1212, 0, 3, 398,
                                                                       408, 738, 753, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1233, 0, 3, 408,
                                                                       418, 753, 768, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1254, 0, 3, 418,
                                                                       428, 768, 783, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1275, 0, 3, 428,
                                                                       438, 783, 798, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1296, 0, 3, 438,
                                                                       448, 798, 813, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1317, 0, 3, 448,
                                                                       458, 813, 828, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1338, 0, 3, 458,
                                                                       468, 828, 843, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1359, 0, 3, 468,
                                                                       478, 843, 858, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1380, 0, 3, 478,
                                                                       488, 858, 873, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1401, 0, 3, 488,
                                                                       498, 873, 888, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1422, 0, 3, 518,
                                                                       528, 903, 918, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1443, 0, 3, 528,
                                                                       538, 918, 933, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1464, 0, 3, 538,
                                                                       548, 933, 948, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1485, 0, 3, 548,
                                                                       558, 948, 963, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1506, 0, 3, 558,
                                                                       568, 963, 978, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1527, 0, 3, 568,
                                                                       578, 978, 993, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1548, 0, 3, 578,
                                                                       588, 993, 1008, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1569, 0, 3, 588,
                                                                       598, 1008, 1023, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1590, 0, 3, 598,
                                                                       608, 1023, 1038, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1611, 0, 3, 608,
                                                                       618, 1038, 1053, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1632, 0, 3, 618,
                                                                       628, 1053, 1068, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1653, 0, 3, 628,
                                                                       638, 1068, 1083, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1674, 0, 3, 638,
                                                                       648, 1083, 1098, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1695, 0, 3, 648,
                                                                       658, 1098, 1113, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1716, 0, 3, 678,
                                                                       693, 1128, 1149, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1744, 0, 3, 693,
                                                                       708, 1149, 1170, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1772, 0, 3, 708,
                                                                       723, 1170, 1191, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1800, 0, 3, 723,
                                                                       738, 1191, 1212, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1828, 0, 3, 738,
                                                                       753, 1212, 1233, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1856, 0, 3, 753,
                                                                       768, 1233, 1254, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1884, 0, 3, 768,
                                                                       783, 1254, 1275, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1912, 0, 3, 783,
                                                                       798, 1275, 1296, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1940, 0, 3, 798,
                                                                       813, 1296, 1317, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1968, 0, 3, 813,
                                                                       828, 1317, 1338, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1996, 0, 3, 828,
                                                                       843, 1338, 1359, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2024, 0, 3, 843,
                                                                       858, 1359, 1380, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2052, 0, 3, 858,
                                                                       873, 1380, 1401, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2080, 0, 3, 903,
                                                                       918, 1422, 1443, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2108, 0, 3, 918,
                                                                       933, 1443, 1464, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2136, 0, 3, 933,
                                                                       948, 1464, 1485, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2164, 0, 3, 948,
                                                                       963, 1485, 1506, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2192, 0, 3, 963,
                                                                       978, 1506, 1527, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2220, 0, 3, 978,
                                                                       993, 1527, 1548, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2248, 0, 3, 993,
                                                                       1008, 1548, 1569, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2276, 0, 3, 1008,
                                                                       1023, 1569, 1590, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2304, 0, 3, 1023,
                                                                       1038, 1590, 1611, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2332, 0, 3, 1038,
                                                                       1053, 1611, 1632, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2360, 0, 3, 1053,
                                                                       1068, 1632, 1653, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2388, 0, 3, 1068,
                                                                       1083, 1653, 1674, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2416, 0, 3, 1083,
                                                                       1098, 1674, 1695, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2444, 0, 3, 1128,
                                                                       1149, 1716, 1744, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2480, 0, 3, 1149,
                                                                       1170, 1744, 1772, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2516, 0, 3, 1170,
                                                                       1191, 1772, 1800, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2552, 0, 3, 1191,
                                                                       1212, 1800, 1828, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2588, 0, 3, 1212,
                                                                       1233, 1828, 1856, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2624, 0, 3, 1233,
                                                                       1254, 1856, 1884, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2660, 0, 3, 1254,
                                                                       1275, 1884, 1912, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2696, 0, 3, 1275,
                                                                       1296, 1912, 1940, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2732, 0, 3, 1296,
                                                                       1317, 1940, 1968, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2768, 0, 3, 1317,
                                                                       1338, 1968, 1996, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2804, 0, 3, 1338,
                                                                       1359, 1996, 2024, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2840, 0, 3, 1359,
                                                                       1380, 2024, 2052, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2876, 0, 3, 1422,
                                                                       1443, 2080, 2108, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2912, 0, 3, 1443,
                                                                       1464, 2108, 2136, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2948, 0, 3, 1464,
                                                                       1485, 2136, 2164, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2984, 0, 3, 1485,
                                                                       1506, 2164, 2192, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 3020, 0, 3, 1506,
                                                                       1527, 2192, 2220, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 3056, 0, 3, 1527,
                                                                       1548, 2220, 2248, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 3092, 0, 3, 1548,
                                                                       1569, 2248, 2276, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 3128, 0, 3, 1569,
                                                                       1590, 2276, 2304, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 3164, 0, 3, 1590,
                                                                       1611, 2304, 2332, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 3200, 0, 3, 1611,
                                                                       1632, 2332, 2360, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 3236, 0, 3, 1632,
                                                                       1653, 2360, 2388, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 3272, 0, 3, 1653,
                                                                       1674, 2388, 2416, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3308, 0, 3, 1716,
                                                                       1744, 2444, 2480, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3353, 0, 3, 1744,
                                                                       1772, 2480, 2516, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3398, 0, 3, 1772,
                                                                       1800, 2516, 2552, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3443, 0, 3, 1800,
                                                                       1828, 2552, 2588, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3488, 0, 3, 1828,
                                                                       1856, 2588, 2624, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3533, 0, 3, 1856,
                                                                       1884, 2624, 2660, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3578, 0, 3, 1884,
                                                                       1912, 2660, 2696, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3623, 0, 3, 1912,
                                                                       1940, 2696, 2732, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3668, 0, 3, 1940,
                                                                       1968, 2732, 2768, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3713, 0, 3, 1968,
                                                                       1996, 2768, 2804, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3758, 0, 3, 1996,
                                                                       2024, 2804, 2840, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3803, 0, 3, 2080,
                                                                       2108, 2876, 2912, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3848, 0, 3, 2108,
                                                                       2136, 2912, 2948, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3893, 0, 3, 2136,
                                                                       2164, 2948, 2984, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3938, 0, 3, 2164,
                                                                       2192, 2984, 3020, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3983, 0, 3, 2192,
                                                                       2220, 3020, 3056, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 4028, 0, 3, 2220,
                                                                       2248, 3056, 3092, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 4073, 0, 3, 2248,
                                                                       2276, 3092, 3128, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 4118, 0, 3, 2276,
                                                                       2304, 3128, 3164, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 4163, 0, 3, 2304,
                                                                       2332, 3164, 3200, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 4208, 0, 3, 2332,
                                                                       2360, 3200, 3236, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 4253, 0, 3, 2360,
                                                                       2388, 3236, 3272, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4298, 0, 3, 2444,
                                                                       2480, 3308, 3353, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4353, 0, 3, 2480,
                                                                       2516, 3353, 3398, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4408, 0, 3, 2516,
                                                                       2552, 3398, 3443, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4463, 0, 3, 2552,
                                                                       2588, 3443, 3488, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4518, 0, 3, 2588,
                                                                       2624, 3488, 3533, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4573, 0, 3, 2624,
                                                                       2660, 3533, 3578, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4628, 0, 3, 2660,
                                                                       2696, 3578, 3623, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4683, 0, 3, 2696,
                                                                       2732, 3623, 3668, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4738, 0, 3, 2732,
                                                                       2768, 3668, 3713, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4793, 0, 3, 2768,
                                                                       2804, 3713, 3758, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4848, 0, 3, 2876,
                                                                       2912, 3803, 3848, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4903, 0, 3, 2912,
                                                                       2948, 3848, 3893, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4958, 0, 3, 2948,
                                                                       2984, 3893, 3938, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 5013, 0, 3, 2984,
                                                                       3020, 3938, 3983, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 5068, 0, 3, 3020,
                                                                       3056, 3983, 4028, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 5123, 0, 3, 3056,
                                                                       3092, 4028, 4073, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 5178, 0, 3, 3092,
                                                                       3128, 4073, 4118, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 5233, 0, 3, 3128,
                                                                       3164, 4118, 4163, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 5288, 0, 3, 3164,
                                                                       3200, 4163, 4208, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 5343, 0, 3, 3200,
                                                                       3236, 4208, 4253, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 5398, 0, 3, 3308,
                                                                       3353, 4298, 4353, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 5464, 0, 3, 3353,
                                                                       3398, 4353, 4408, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 5530, 0, 3, 3398,
                                                                       3443, 4408, 4463, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 5596, 0, 3, 3443,
                                                                       3488, 4463, 4518, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 5662, 0, 3, 3488,
                                                                       3533, 4518, 4573, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 5728, 0, 3, 3533,
                                                                       3578, 4573, 4628, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 5794, 0, 3, 3578,
                                                                       3623, 4628, 4683, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 5860, 0, 3, 3623,
                                                                       3668, 4683, 4738, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 5926, 0, 3, 3668,
                                                                       3713, 4738, 4793, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 5992, 0, 3, 3803,
                                                                       3848, 4848, 4903, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 6058, 0, 3, 3848,
                                                                       3893, 4903, 4958, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 6124, 0, 3, 3893,
                                                                       3938, 4958, 5013, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 6190, 0, 3, 3938,
                                                                       3983, 5013, 5068, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 6256, 0, 3, 3983,
                                                                       4028, 5068, 5123, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 6322, 0, 3, 4028,
                                                                       4073, 5123, 5178, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 6388, 0, 3, 4073,
                                                                       4118, 5178, 5233, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 6454, 0, 3, 4118,
                                                                       4163, 5233, 5288, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 6520, 0, 3, 4163,
                                                                       4208, 5288, 5343, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 6586, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 6589, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 6592, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 6595, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 6598, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 6601, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 6604, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 6607, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 6610, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 6613, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 6616, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 6619, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 6622, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 6625, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 6628, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 6631, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 6634, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 6637, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 6640, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 6643, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 6646, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 6649, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 6652, 3, 34,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 6655, 3, 35,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 6658, 3, 36,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 6661, 3, 37,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 6664, 3, 38,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 6667, 3, 39,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 6670, 3, 40,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 6673, 3, 41,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 6676, 3, 42,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 6679, 3, 43,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 6682, 3, 44,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 6685, 3, 45,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 6688, 3, 9, 52,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 6697, 3, 10, 55,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 6706, 3, 11, 58,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 6715, 3, 12, 61,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 6724, 3, 13, 64,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 6733, 3, 14, 67,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 6742, 3, 15, 70,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 6751, 3, 16, 73,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 6760, 3, 17, 76,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 6769, 3, 18, 79,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 6778, 3, 19, 82,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 6787, 3, 20, 85,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 6796, 3, 21, 88,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 6805, 3, 22, 91,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 6814, 3, 23, 94,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 6823, 3, 24, 97,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 6832, 3, 29, 106,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 6841, 3, 30, 109,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 6850, 3, 31, 112,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 6859, 3, 32, 115,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 6868, 3, 33, 118,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 6877, 3, 34, 121,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 6886, 3, 35, 124,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 6895, 3, 36, 127,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 6904, 3, 37, 130,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 6913, 3, 38, 133,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 6922, 3, 39, 136,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 6931, 3, 40, 139,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 6940, 3, 41, 142,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 6949, 3, 42, 145,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 6958, 3, 43, 148,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 6967, 3, 44, 151,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6976, 3, 52, 166,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6994, 3, 55, 172,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7012, 3, 58, 178,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7030, 3, 61, 184,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7048, 3, 64, 190,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7066, 3, 67, 196,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7084, 3, 70, 202,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7102, 3, 73, 208,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7120, 3, 76, 214,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7138, 3, 79, 220,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7156, 3, 82, 226,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7174, 3, 85, 232,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7192, 3, 88, 238,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7210, 3, 91, 244,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7228, 3, 94, 250,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7246, 3, 106, 268,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7264, 3, 109, 274,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7282, 3, 112, 280,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7300, 3, 115, 286,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7318, 3, 118, 292,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7336, 3, 121, 298,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7354, 3, 124, 304,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7372, 3, 127, 310,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7390, 3, 130, 316,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7408, 3, 133, 322,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7426, 3, 136, 328,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7444, 3, 139, 334,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7462, 3, 142, 340,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7480, 3, 145, 346,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7498, 3, 148, 352,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 7516, 3, 166, 378,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 7546, 3, 172, 388,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 7576, 3, 178, 398,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 7606, 3, 184, 408,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 7636, 3, 190, 418,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 7666, 3, 196, 428,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 7696, 3, 202, 438,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 7726, 3, 208, 448,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 7756, 3, 214, 458,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 7786, 3, 220, 468,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 7816, 3, 226, 478,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 7846, 3, 232, 488,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 7876, 3, 238, 498,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 7906, 3, 244, 508,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 7936, 3, 268, 538,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 7966, 3, 274, 548,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 7996, 3, 280, 558,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 8026, 3, 286, 568,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 8056, 3, 292, 578,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 8086, 3, 298, 588,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 8116, 3, 304, 598,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 8146, 3, 310, 608,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 8176, 3, 316, 618,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 8206, 3, 322, 628,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 8236, 3, 328, 638,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 8266, 3, 334, 648,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 8296, 3, 340, 658,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 8326, 3, 346, 668,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 8356, 3, 378, 708,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 8401, 3, 388, 723,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 8446, 3, 398, 738,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 8491, 3, 408, 753,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 8536, 3, 418, 768,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 8581, 3, 428, 783,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 8626, 3, 438, 798,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 8671, 3, 448, 813,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 8716, 3, 458, 828,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 8761, 3, 468, 843,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 8806, 3, 478, 858,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 8851, 3, 488, 873,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 8896, 3, 498, 888,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 8941, 3, 538, 933,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 8986, 3, 548, 948,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 9031, 3, 558, 963,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 9076, 3, 568, 978,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 9121, 3, 578, 993,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 9166, 3, 588,
                                                                       1008, ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 9211, 3, 598,
                                                                       1023, ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 9256, 3, 608,
                                                                       1038, ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 9301, 3, 618,
                                                                       1053, ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 9346, 3, 628,
                                                                       1068, ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 9391, 3, 638,
                                                                       1083, ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 9436, 3, 648,
                                                                       1098, ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 9481, 3, 658,
                                                                       1113, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 9526, 3, 708,
                                                                       1170, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 9589, 3, 723,
                                                                       1191, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 9652, 3, 738,
                                                                       1212, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 9715, 3, 753,
                                                                       1233, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 9778, 3, 768,
                                                                       1254, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 9841, 3, 783,
                                                                       1275, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 9904, 3, 798,
                                                                       1296, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 9967, 3, 813,
                                                                       1317, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 10030, 3, 828,
                                                                       1338, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 10093, 3, 843,
                                                                       1359, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 10156, 3, 858,
                                                                       1380, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 10219, 3, 873,
                                                                       1401, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 10282, 3, 933,
                                                                       1464, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 10345, 3, 948,
                                                                       1485, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 10408, 3, 963,
                                                                       1506, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 10471, 3, 978,
                                                                       1527, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 10534, 3, 993,
                                                                       1548, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 10597, 3, 1008,
                                                                       1569, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 10660, 3, 1023,
                                                                       1590, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 10723, 3, 1038,
                                                                       1611, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 10786, 3, 1053,
                                                                       1632, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 10849, 3, 1068,
                                                                       1653, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 10912, 3, 1083,
                                                                       1674, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 10975, 3, 1098,
                                                                       1695, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 11038, 3, 1170,
                                                                       1772, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 11122, 3, 1191,
                                                                       1800, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 11206, 3, 1212,
                                                                       1828, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 11290, 3, 1233,
                                                                       1856, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 11374, 3, 1254,
                                                                       1884, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 11458, 3, 1275,
                                                                       1912, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 11542, 3, 1296,
                                                                       1940, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 11626, 3, 1317,
                                                                       1968, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 11710, 3, 1338,
                                                                       1996, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 11794, 3, 1359,
                                                                       2024, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 11878, 3, 1380,
                                                                       2052, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 11962, 3, 1464,
                                                                       2136, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 12046, 3, 1485,
                                                                       2164, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 12130, 3, 1506,
                                                                       2192, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 12214, 3, 1527,
                                                                       2220, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 12298, 3, 1548,
                                                                       2248, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 12382, 3, 1569,
                                                                       2276, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 12466, 3, 1590,
                                                                       2304, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 12550, 3, 1611,
                                                                       2332, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 12634, 3, 1632,
                                                                       2360, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 12718, 3, 1653,
                                                                       2388, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 12802, 3, 1674,
                                                                       2416, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 12886, 3, 1772,
                                                                       2516, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 12994, 3, 1800,
                                                                       2552, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 13102, 3, 1828,
                                                                       2588, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 13210, 3, 1856,
                                                                       2624, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 13318, 3, 1884,
                                                                       2660, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 13426, 3, 1912,
                                                                       2696, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 13534, 3, 1940,
                                                                       2732, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 13642, 3, 1968,
                                                                       2768, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 13750, 3, 1996,
                                                                       2804, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 13858, 3, 2024,
                                                                       2840, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 13966, 3, 2136,
                                                                       2948, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 14074, 3, 2164,
                                                                       2984, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 14182, 3, 2192,
                                                                       3020, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 14290, 3, 2220,
                                                                       3056, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 14398, 3, 2248,
                                                                       3092, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 14506, 3, 2276,
                                                                       3128, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 14614, 3, 2304,
                                                                       3164, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 14722, 3, 2332,
                                                                       3200, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 14830, 3, 2360,
                                                                       3236, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 14938, 3, 2388,
                                                                       3272, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 15046, 3, 2516,
                                                                       3398, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 15181, 3, 2552,
                                                                       3443, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 15316, 3, 2588,
                                                                       3488, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 15451, 3, 2624,
                                                                       3533, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 15586, 3, 2660,
                                                                       3578, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 15721, 3, 2696,
                                                                       3623, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 15856, 3, 2732,
                                                                       3668, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 15991, 3, 2768,
                                                                       3713, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 16126, 3, 2804,
                                                                       3758, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 16261, 3, 2948,
                                                                       3893, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 16396, 3, 2984,
                                                                       3938, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 16531, 3, 3020,
                                                                       3983, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 16666, 3, 3056,
                                                                       4028, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 16801, 3, 3092,
                                                                       4073, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 16936, 3, 3128,
                                                                       4118, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 17071, 3, 3164,
                                                                       4163, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 17206, 3, 3200,
                                                                       4208, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 17341, 3, 3236,
                                                                       4253, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 17476, 3, 3398,
                                                                       4408, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 17641, 3, 3443,
                                                                       4463, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 17806, 3, 3488,
                                                                       4518, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 17971, 3, 3533,
                                                                       4573, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 18136, 3, 3578,
                                                                       4628, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 18301, 3, 3623,
                                                                       4683, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 18466, 3, 3668,
                                                                       4738, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 18631, 3, 3713,
                                                                       4793, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 18796, 3, 3893,
                                                                       4958, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 18961, 3, 3938,
                                                                       5013, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 19126, 3, 3983,
                                                                       5068, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 19291, 3, 4028,
                                                                       5123, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 19456, 3, 4073,
                                                                       5178, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 19621, 3, 4118,
                                                                       5233, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 19786, 3, 4163,
                                                                       5288, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 19951, 3, 4208,
                                                                       5343, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 20116, 3, 4408,
                                                                       5530, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 20314, 3, 4463,
                                                                       5596, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 20512, 3, 4518,
                                                                       5662, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 20710, 3, 4573,
                                                                       5728, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 20908, 3, 4628,
                                                                       5794, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 21106, 3, 4683,
                                                                       5860, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 21304, 3, 4738,
                                                                       5926, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 21502, 3, 4958,
                                                                       6124, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 21700, 3, 5013,
                                                                       6190, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 21898, 3, 5068,
                                                                       6256, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 22096, 3, 5123,
                                                                       6322, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 22294, 3, 5178,
                                                                       6388, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 22492, 3, 5233,
                                                                       6454, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 22690, 3, 5288,
                                                                       6520, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22888, 3, 7, 8,
                                                                       6586, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22894, 3, 8, 9,
                                                                       6589, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22900, 3, 9, 10,
                                                                       6592, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22906, 3, 10, 11,
                                                                       6595, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22912, 3, 11, 12,
                                                                       6598, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22918, 3, 12, 13,
                                                                       6601, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22924, 3, 13, 14,
                                                                       6604, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22930, 3, 14, 15,
                                                                       6607, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22936, 3, 15, 16,
                                                                       6610, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22942, 3, 16, 17,
                                                                       6613, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22948, 3, 17, 18,
                                                                       6616, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22954, 3, 18, 19,
                                                                       6619, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22960, 3, 19, 20,
                                                                       6622, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22966, 3, 20, 21,
                                                                       6625, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22972, 3, 21, 22,
                                                                       6628, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22978, 3, 22, 23,
                                                                       6631, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22984, 3, 23, 24,
                                                                       6634, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22990, 3, 27, 28,
                                                                       6637, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22996, 3, 28, 29,
                                                                       6640, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23002, 3, 29, 30,
                                                                       6643, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23008, 3, 30, 31,
                                                                       6646, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23014, 3, 31, 32,
                                                                       6649, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23020, 3, 32, 33,
                                                                       6652, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23026, 3, 33, 34,
                                                                       6655, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23032, 3, 34, 35,
                                                                       6658, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23038, 3, 35, 36,
                                                                       6661, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23044, 3, 36, 37,
                                                                       6664, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23050, 3, 37, 38,
                                                                       6667, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23056, 3, 38, 39,
                                                                       6670, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23062, 3, 39, 40,
                                                                       6673, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23068, 3, 40, 41,
                                                                       6676, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23074, 3, 41, 42,
                                                                       6679, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23080, 3, 42, 43,
                                                                       6682, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23086, 3, 43, 44,
                                                                       6685, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 23092, 0, 3,
                                                                       22888, 6586, 22894, 6688,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 23110, 0, 3,
                                                                       22894, 6589, 22900, 6697,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 23128, 0, 3,
                                                                       22900, 6592, 22906, 6706,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 23146, 0, 3,
                                                                       22906, 6595, 22912, 6715,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 23164, 0, 3,
                                                                       22912, 6598, 22918, 6724,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 23182, 0, 3,
                                                                       22918, 6601, 22924, 6733,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 23200, 0, 3,
                                                                       22924, 6604, 22930, 6742,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 23218, 0, 3,
                                                                       22930, 6607, 22936, 6751,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 23236, 0, 3,
                                                                       22936, 6610, 22942, 6760,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 23254, 0, 3,
                                                                       22942, 6613, 22948, 6769,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 23272, 0, 3,
                                                                       22948, 6616, 22954, 6778,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 23290, 0, 3,
                                                                       22954, 6619, 22960, 6787,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 23308, 0, 3,
                                                                       22960, 6622, 22966, 6796,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 23326, 0, 3,
                                                                       22966, 6625, 22972, 6805,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 23344, 0, 3,
                                                                       22972, 6628, 22978, 6814,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 23362, 0, 3,
                                                                       22978, 6631, 22984, 6823,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 23380, 0, 3,
                                                                       22990, 6637, 22996, 6832,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 23398, 0, 3,
                                                                       22996, 6640, 23002, 6841,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 23416, 0, 3,
                                                                       23002, 6643, 23008, 6850,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 23434, 0, 3,
                                                                       23008, 6646, 23014, 6859,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 23452, 0, 3,
                                                                       23014, 6649, 23020, 6868,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 23470, 0, 3,
                                                                       23020, 6652, 23026, 6877,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 23488, 0, 3,
                                                                       23026, 6655, 23032, 6886,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 23506, 0, 3,
                                                                       23032, 6658, 23038, 6895,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 23524, 0, 3,
                                                                       23038, 6661, 23044, 6904,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 23542, 0, 3,
                                                                       23044, 6664, 23050, 6913,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 23560, 0, 3,
                                                                       23050, 6667, 23056, 6922,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 23578, 0, 3,
                                                                       23056, 6670, 23062, 6931,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 23596, 0, 3,
                                                                       23062, 6673, 23068, 6940,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 23614, 0, 3,
                                                                       23068, 6676, 23074, 6949,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 23632, 0, 3,
                                                                       23074, 6679, 23080, 6958,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 23650, 0, 3,
                                                                       23080, 6682, 23086, 6967,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 23668, 0, 3,
                                                                       23092, 6688, 23110, 154,
                                                                       160, 6976, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 23704, 0, 3,
                                                                       23110, 6697, 23128, 160,
                                                                       166, 6994, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 23740, 0, 3,
                                                                       23128, 6706, 23146, 166,
                                                                       172, 7012, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 23776, 0, 3,
                                                                       23146, 6715, 23164, 172,
                                                                       178, 7030, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 23812, 0, 3,
                                                                       23164, 6724, 23182, 178,
                                                                       184, 7048, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 23848, 0, 3,
                                                                       23182, 6733, 23200, 184,
                                                                       190, 7066, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 23884, 0, 3,
                                                                       23200, 6742, 23218, 190,
                                                                       196, 7084, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 23920, 0, 3,
                                                                       23218, 6751, 23236, 196,
                                                                       202, 7102, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 23956, 0, 3,
                                                                       23236, 6760, 23254, 202,
                                                                       208, 7120, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 23992, 0, 3,
                                                                       23254, 6769, 23272, 208,
                                                                       214, 7138, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 24028, 0, 3,
                                                                       23272, 6778, 23290, 214,
                                                                       220, 7156, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 24064, 0, 3,
                                                                       23290, 6787, 23308, 220,
                                                                       226, 7174, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 24100, 0, 3,
                                                                       23308, 6796, 23326, 226,
                                                                       232, 7192, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 24136, 0, 3,
                                                                       23326, 6805, 23344, 232,
                                                                       238, 7210, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 24172, 0, 3,
                                                                       23344, 6814, 23362, 238,
                                                                       244, 7228, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 24208, 0, 3,
                                                                       23380, 6832, 23398, 256,
                                                                       262, 7246, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 24244, 0, 3,
                                                                       23398, 6841, 23416, 262,
                                                                       268, 7264, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 24280, 0, 3,
                                                                       23416, 6850, 23434, 268,
                                                                       274, 7282, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 24316, 0, 3,
                                                                       23434, 6859, 23452, 274,
                                                                       280, 7300, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 24352, 0, 3,
                                                                       23452, 6868, 23470, 280,
                                                                       286, 7318, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 24388, 0, 3,
                                                                       23470, 6877, 23488, 286,
                                                                       292, 7336, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 24424, 0, 3,
                                                                       23488, 6886, 23506, 292,
                                                                       298, 7354, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 24460, 0, 3,
                                                                       23506, 6895, 23524, 298,
                                                                       304, 7372, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 24496, 0, 3,
                                                                       23524, 6904, 23542, 304,
                                                                       310, 7390, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 24532, 0, 3,
                                                                       23542, 6913, 23560, 310,
                                                                       316, 7408, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 24568, 0, 3,
                                                                       23560, 6922, 23578, 316,
                                                                       322, 7426, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 24604, 0, 3,
                                                                       23578, 6931, 23596, 322,
                                                                       328, 7444, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 24640, 0, 3,
                                                                       23596, 6940, 23614, 328,
                                                                       334, 7462, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 24676, 0, 3,
                                                                       23614, 6949, 23632, 334,
                                                                       340, 7480, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 24712, 0, 3,
                                                                       23632, 6958, 23650, 340,
                                                                       346, 7498, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 24748, 0, 3,
                                                                       23668, 6976, 23704, 358,
                                                                       368, 7516, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 24808, 0, 3,
                                                                       23704, 6994, 23740, 368,
                                                                       378, 7546, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 24868, 0, 3,
                                                                       23740, 7012, 23776, 378,
                                                                       388, 7576, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 24928, 0, 3,
                                                                       23776, 7030, 23812, 388,
                                                                       398, 7606, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 24988, 0, 3,
                                                                       23812, 7048, 23848, 398,
                                                                       408, 7636, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 25048, 0, 3,
                                                                       23848, 7066, 23884, 408,
                                                                       418, 7666, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 25108, 0, 3,
                                                                       23884, 7084, 23920, 418,
                                                                       428, 7696, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 25168, 0, 3,
                                                                       23920, 7102, 23956, 428,
                                                                       438, 7726, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 25228, 0, 3,
                                                                       23956, 7120, 23992, 438,
                                                                       448, 7756, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 25288, 0, 3,
                                                                       23992, 7138, 24028, 448,
                                                                       458, 7786, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 25348, 0, 3,
                                                                       24028, 7156, 24064, 458,
                                                                       468, 7816, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 25408, 0, 3,
                                                                       24064, 7174, 24100, 468,
                                                                       478, 7846, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 25468, 0, 3,
                                                                       24100, 7192, 24136, 478,
                                                                       488, 7876, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 25528, 0, 3,
                                                                       24136, 7210, 24172, 488,
                                                                       498, 7906, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 25588, 0, 3,
                                                                       24208, 7246, 24244, 518,
                                                                       528, 7936, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 25648, 0, 3,
                                                                       24244, 7264, 24280, 528,
                                                                       538, 7966, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 25708, 0, 3,
                                                                       24280, 7282, 24316, 538,
                                                                       548, 7996, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 25768, 0, 3,
                                                                       24316, 7300, 24352, 548,
                                                                       558, 8026, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 25828, 0, 3,
                                                                       24352, 7318, 24388, 558,
                                                                       568, 8056, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 25888, 0, 3,
                                                                       24388, 7336, 24424, 568,
                                                                       578, 8086, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 25948, 0, 3,
                                                                       24424, 7354, 24460, 578,
                                                                       588, 8116, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 26008, 0, 3,
                                                                       24460, 7372, 24496, 588,
                                                                       598, 8146, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 26068, 0, 3,
                                                                       24496, 7390, 24532, 598,
                                                                       608, 8176, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 26128, 0, 3,
                                                                       24532, 7408, 24568, 608,
                                                                       618, 8206, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 26188, 0, 3,
                                                                       24568, 7426, 24604, 618,
                                                                       628, 8236, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 26248, 0, 3,
                                                                       24604, 7444, 24640, 628,
                                                                       638, 8266, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 26308, 0, 3,
                                                                       24640, 7462, 24676, 638,
                                                                       648, 8296, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 26368, 0, 3,
                                                                       24676, 7480, 24712, 648,
                                                                       658, 8326, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 26428, 0, 3,
                                                                       24748, 7516, 24808, 678,
                                                                       693, 8356, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 26518, 0, 3,
                                                                       24808, 7546, 24868, 693,
                                                                       708, 8401, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 26608, 0, 3,
                                                                       24868, 7576, 24928, 708,
                                                                       723, 8446, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 26698, 0, 3,
                                                                       24928, 7606, 24988, 723,
                                                                       738, 8491, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 26788, 0, 3,
                                                                       24988, 7636, 25048, 738,
                                                                       753, 8536, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 26878, 0, 3,
                                                                       25048, 7666, 25108, 753,
                                                                       768, 8581, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 26968, 0, 3,
                                                                       25108, 7696, 25168, 768,
                                                                       783, 8626, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 27058, 0, 3,
                                                                       25168, 7726, 25228, 783,
                                                                       798, 8671, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 27148, 0, 3,
                                                                       25228, 7756, 25288, 798,
                                                                       813, 8716, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 27238, 0, 3,
                                                                       25288, 7786, 25348, 813,
                                                                       828, 8761, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 27328, 0, 3,
                                                                       25348, 7816, 25408, 828,
                                                                       843, 8806, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 27418, 0, 3,
                                                                       25408, 7846, 25468, 843,
                                                                       858, 8851, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 27508, 0, 3,
                                                                       25468, 7876, 25528, 858,
                                                                       873, 8896, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 27598, 0, 3,
                                                                       25588, 7936, 25648, 903,
                                                                       918, 8941, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 27688, 0, 3,
                                                                       25648, 7966, 25708, 918,
                                                                       933, 8986, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 27778, 0, 3,
                                                                       25708, 7996, 25768, 933,
                                                                       948, 9031, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 27868, 0, 3,
                                                                       25768, 8026, 25828, 948,
                                                                       963, 9076, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 27958, 0, 3,
                                                                       25828, 8056, 25888, 963,
                                                                       978, 9121, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 28048, 0, 3,
                                                                       25888, 8086, 25948, 978,
                                                                       993, 9166, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 28138, 0, 3,
                                                                       25948, 8116, 26008, 993,
                                                                       1008, 9211, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 28228, 0, 3,
                                                                       26008, 8146, 26068, 1008,
                                                                       1023, 9256, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 28318, 0, 3,
                                                                       26068, 8176, 26128, 1023,
                                                                       1038, 9301, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 28408, 0, 3,
                                                                       26128, 8206, 26188, 1038,
                                                                       1053, 9346, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 28498, 0, 3,
                                                                       26188, 8236, 26248, 1053,
                                                                       1068, 9391, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 28588, 0, 3,
                                                                       26248, 8266, 26308, 1068,
                                                                       1083, 9436, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 28678, 0, 3,
                                                                       26308, 8296, 26368, 1083,
                                                                       1098, 9481, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 28768, 0, 3,
                                                                       26428, 8356, 26518, 1128,
                                                                       1149, 9526, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 28894, 0, 3,
                                                                       26518, 8401, 26608, 1149,
                                                                       1170, 9589, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 29020, 0, 3,
                                                                       26608, 8446, 26698, 1170,
                                                                       1191, 9652, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 29146, 0, 3,
                                                                       26698, 8491, 26788, 1191,
                                                                       1212, 9715, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 29272, 0, 3,
                                                                       26788, 8536, 26878, 1212,
                                                                       1233, 9778, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 29398, 0, 3,
                                                                       26878, 8581, 26968, 1233,
                                                                       1254, 9841, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 29524, 0, 3,
                                                                       26968, 8626, 27058, 1254,
                                                                       1275, 9904, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 29650, 0, 3,
                                                                       27058, 8671, 27148, 1275,
                                                                       1296, 9967, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 29776, 0, 3,
                                                                       27148, 8716, 27238, 1296,
                                                                       1317, 10030, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 29902, 0, 3,
                                                                       27238, 8761, 27328, 1317,
                                                                       1338, 10093, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 30028, 0, 3,
                                                                       27328, 8806, 27418, 1338,
                                                                       1359, 10156, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 30154, 0, 3,
                                                                       27418, 8851, 27508, 1359,
                                                                       1380, 10219, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 30280, 0, 3,
                                                                       27598, 8941, 27688, 1422,
                                                                       1443, 10282, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 30406, 0, 3,
                                                                       27688, 8986, 27778, 1443,
                                                                       1464, 10345, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 30532, 0, 3,
                                                                       27778, 9031, 27868, 1464,
                                                                       1485, 10408, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 30658, 0, 3,
                                                                       27868, 9076, 27958, 1485,
                                                                       1506, 10471, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 30784, 0, 3,
                                                                       27958, 9121, 28048, 1506,
                                                                       1527, 10534, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 30910, 0, 3,
                                                                       28048, 9166, 28138, 1527,
                                                                       1548, 10597, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 31036, 0, 3,
                                                                       28138, 9211, 28228, 1548,
                                                                       1569, 10660, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 31162, 0, 3,
                                                                       28228, 9256, 28318, 1569,
                                                                       1590, 10723, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 31288, 0, 3,
                                                                       28318, 9301, 28408, 1590,
                                                                       1611, 10786, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 31414, 0, 3,
                                                                       28408, 9346, 28498, 1611,
                                                                       1632, 10849, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 31540, 0, 3,
                                                                       28498, 9391, 28588, 1632,
                                                                       1653, 10912, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 31666, 0, 3,
                                                                       28588, 9436, 28678, 1653,
                                                                       1674, 10975, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 31792, 0, 3,
                                                                       28768, 9526, 28894, 1716,
                                                                       1744, 11038, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 31960, 0, 3,
                                                                       28894, 9589, 29020, 1744,
                                                                       1772, 11122, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 32128, 0, 3,
                                                                       29020, 9652, 29146, 1772,
                                                                       1800, 11206, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 32296, 0, 3,
                                                                       29146, 9715, 29272, 1800,
                                                                       1828, 11290, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 32464, 0, 3,
                                                                       29272, 9778, 29398, 1828,
                                                                       1856, 11374, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 32632, 0, 3,
                                                                       29398, 9841, 29524, 1856,
                                                                       1884, 11458, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 32800, 0, 3,
                                                                       29524, 9904, 29650, 1884,
                                                                       1912, 11542, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 32968, 0, 3,
                                                                       29650, 9967, 29776, 1912,
                                                                       1940, 11626, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 33136, 0, 3,
                                                                       29776, 10030, 29902, 1940,
                                                                       1968, 11710, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 33304, 0, 3,
                                                                       29902, 10093, 30028, 1968,
                                                                       1996, 11794, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 33472, 0, 3,
                                                                       30028, 10156, 30154, 1996,
                                                                       2024, 11878, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 33640, 0, 3,
                                                                       30280, 10282, 30406, 2080,
                                                                       2108, 11962, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 33808, 0, 3,
                                                                       30406, 10345, 30532, 2108,
                                                                       2136, 12046, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 33976, 0, 3,
                                                                       30532, 10408, 30658, 2136,
                                                                       2164, 12130, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 34144, 0, 3,
                                                                       30658, 10471, 30784, 2164,
                                                                       2192, 12214, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 34312, 0, 3,
                                                                       30784, 10534, 30910, 2192,
                                                                       2220, 12298, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 34480, 0, 3,
                                                                       30910, 10597, 31036, 2220,
                                                                       2248, 12382, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 34648, 0, 3,
                                                                       31036, 10660, 31162, 2248,
                                                                       2276, 12466, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 34816, 0, 3,
                                                                       31162, 10723, 31288, 2276,
                                                                       2304, 12550, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 34984, 0, 3,
                                                                       31288, 10786, 31414, 2304,
                                                                       2332, 12634, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 35152, 0, 3,
                                                                       31414, 10849, 31540, 2332,
                                                                       2360, 12718, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 35320, 0, 3,
                                                                       31540, 10912, 31666, 2360,
                                                                       2388, 12802, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 35488, 0, 3,
                                                                       31792, 11038, 31960, 2444,
                                                                       2480, 12886, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 35704, 0, 3,
                                                                       31960, 11122, 32128, 2480,
                                                                       2516, 12994, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 35920, 0, 3,
                                                                       32128, 11206, 32296, 2516,
                                                                       2552, 13102, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 36136, 0, 3,
                                                                       32296, 11290, 32464, 2552,
                                                                       2588, 13210, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 36352, 0, 3,
                                                                       32464, 11374, 32632, 2588,
                                                                       2624, 13318, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 36568, 0, 3,
                                                                       32632, 11458, 32800, 2624,
                                                                       2660, 13426, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 36784, 0, 3,
                                                                       32800, 11542, 32968, 2660,
                                                                       2696, 13534, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 37000, 0, 3,
                                                                       32968, 11626, 33136, 2696,
                                                                       2732, 13642, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 37216, 0, 3,
                                                                       33136, 11710, 33304, 2732,
                                                                       2768, 13750, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 37432, 0, 3,
                                                                       33304, 11794, 33472, 2768,
                                                                       2804, 13858, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 37648, 0, 3,
                                                                       33640, 11962, 33808, 2876,
                                                                       2912, 13966, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 37864, 0, 3,
                                                                       33808, 12046, 33976, 2912,
                                                                       2948, 14074, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 38080, 0, 3,
                                                                       33976, 12130, 34144, 2948,
                                                                       2984, 14182, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 38296, 0, 3,
                                                                       34144, 12214, 34312, 2984,
                                                                       3020, 14290, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 38512, 0, 3,
                                                                       34312, 12298, 34480, 3020,
                                                                       3056, 14398, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 38728, 0, 3,
                                                                       34480, 12382, 34648, 3056,
                                                                       3092, 14506, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 38944, 0, 3,
                                                                       34648, 12466, 34816, 3092,
                                                                       3128, 14614, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 39160, 0, 3,
                                                                       34816, 12550, 34984, 3128,
                                                                       3164, 14722, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 39376, 0, 3,
                                                                       34984, 12634, 35152, 3164,
                                                                       3200, 14830, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 39592, 0, 3,
                                                                       35152, 12718, 35320, 3200,
                                                                       3236, 14938, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 39808, 0, 3,
                                                                       35488, 12886, 35704, 3308,
                                                                       3353, 15046, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 40078, 0, 3,
                                                                       35704, 12994, 35920, 3353,
                                                                       3398, 15181, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 40348, 0, 3,
                                                                       35920, 13102, 36136, 3398,
                                                                       3443, 15316, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 40618, 0, 3,
                                                                       36136, 13210, 36352, 3443,
                                                                       3488, 15451, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 40888, 0, 3,
                                                                       36352, 13318, 36568, 3488,
                                                                       3533, 15586, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 41158, 0, 3,
                                                                       36568, 13426, 36784, 3533,
                                                                       3578, 15721, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 41428, 0, 3,
                                                                       36784, 13534, 37000, 3578,
                                                                       3623, 15856, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 41698, 0, 3,
                                                                       37000, 13642, 37216, 3623,
                                                                       3668, 15991, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 41968, 0, 3,
                                                                       37216, 13750, 37432, 3668,
                                                                       3713, 16126, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 42238, 0, 3,
                                                                       37648, 13966, 37864, 3803,
                                                                       3848, 16261, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 42508, 0, 3,
                                                                       37864, 14074, 38080, 3848,
                                                                       3893, 16396, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 42778, 0, 3,
                                                                       38080, 14182, 38296, 3893,
                                                                       3938, 16531, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 43048, 0, 3,
                                                                       38296, 14290, 38512, 3938,
                                                                       3983, 16666, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 43318, 0, 3,
                                                                       38512, 14398, 38728, 3983,
                                                                       4028, 16801, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 43588, 0, 3,
                                                                       38728, 14506, 38944, 4028,
                                                                       4073, 16936, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 43858, 0, 3,
                                                                       38944, 14614, 39160, 4073,
                                                                       4118, 17071, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 44128, 0, 3,
                                                                       39160, 14722, 39376, 4118,
                                                                       4163, 17206, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 44398, 0, 3,
                                                                       39376, 14830, 39592, 4163,
                                                                       4208, 17341, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 44668, 0, 3,
                                                                       39808, 15046, 40078, 4298,
                                                                       4353, 17476, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 44998, 0, 3,
                                                                       40078, 15181, 40348, 4353,
                                                                       4408, 17641, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 45328, 0, 3,
                                                                       40348, 15316, 40618, 4408,
                                                                       4463, 17806, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 45658, 0, 3,
                                                                       40618, 15451, 40888, 4463,
                                                                       4518, 17971, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 45988, 0, 3,
                                                                       40888, 15586, 41158, 4518,
                                                                       4573, 18136, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 46318, 0, 3,
                                                                       41158, 15721, 41428, 4573,
                                                                       4628, 18301, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 46648, 0, 3,
                                                                       41428, 15856, 41698, 4628,
                                                                       4683, 18466, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 46978, 0, 3,
                                                                       41698, 15991, 41968, 4683,
                                                                       4738, 18631, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 47308, 0, 3,
                                                                       42238, 16261, 42508, 4848,
                                                                       4903, 18796, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 47638, 0, 3,
                                                                       42508, 16396, 42778, 4903,
                                                                       4958, 18961, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 47968, 0, 3,
                                                                       42778, 16531, 43048, 4958,
                                                                       5013, 19126, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 48298, 0, 3,
                                                                       43048, 16666, 43318, 5013,
                                                                       5068, 19291, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 48628, 0, 3,
                                                                       43318, 16801, 43588, 5068,
                                                                       5123, 19456, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 48958, 0, 3,
                                                                       43588, 16936, 43858, 5123,
                                                                       5178, 19621, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 49288, 0, 3,
                                                                       43858, 17071, 44128, 5178,
                                                                       5233, 19786, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 49618, 0, 3,
                                                                       44128, 17206, 44398, 5233,
                                                                       5288, 19951, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 49948, 0, 3,
                                                                       44668, 17476, 44998, 5398,
                                                                       5464, 20116, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 50344, 0, 3,
                                                                       44998, 17641, 45328, 5464,
                                                                       5530, 20314, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 50740, 0, 3,
                                                                       45328, 17806, 45658, 5530,
                                                                       5596, 20512, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 51136, 0, 3,
                                                                       45658, 17971, 45988, 5596,
                                                                       5662, 20710, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 51532, 0, 3,
                                                                       45988, 18136, 46318, 5662,
                                                                       5728, 20908, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 51928, 0, 3,
                                                                       46318, 18301, 46648, 5728,
                                                                       5794, 21106, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 52324, 0, 3,
                                                                       46648, 18466, 46978, 5794,
                                                                       5860, 21304, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 52720, 0, 3,
                                                                       47308, 18796, 47638, 5992,
                                                                       6058, 21502, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 53116, 0, 3,
                                                                       47638, 18961, 47968, 6058,
                                                                       6124, 21700, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 53512, 0, 3,
                                                                       47968, 19126, 48298, 6124,
                                                                       6190, 21898, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 53908, 0, 3,
                                                                       48298, 19291, 48628, 6190,
                                                                       6256, 22096, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 54304, 0, 3,
                                                                       48628, 19456, 48958, 6256,
                                                                       6322, 22294, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 54700, 0, 3,
                                                                       48958, 19621, 49288, 6322,
                                                                       6388, 22492, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 55096, 0, 3,
                                                                       49288, 19786, 49618, 6388,
                                                                       6454, 22690, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 55492, 3, 6586,
                                                                       6589, 22900, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 55502, 3, 6589,
                                                                       6592, 22906, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 55512, 3, 6592,
                                                                       6595, 22912, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 55522, 3, 6595,
                                                                       6598, 22918, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 55532, 3, 6598,
                                                                       6601, 22924, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 55542, 3, 6601,
                                                                       6604, 22930, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 55552, 3, 6604,
                                                                       6607, 22936, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 55562, 3, 6607,
                                                                       6610, 22942, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 55572, 3, 6610,
                                                                       6613, 22948, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 55582, 3, 6613,
                                                                       6616, 22954, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 55592, 3, 6616,
                                                                       6619, 22960, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 55602, 3, 6619,
                                                                       6622, 22966, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 55612, 3, 6622,
                                                                       6625, 22972, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 55622, 3, 6625,
                                                                       6628, 22978, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 55632, 3, 6628,
                                                                       6631, 22984, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 55642, 3, 6637,
                                                                       6640, 23002, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 55652, 3, 6640,
                                                                       6643, 23008, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 55662, 3, 6643,
                                                                       6646, 23014, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 55672, 3, 6646,
                                                                       6649, 23020, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 55682, 3, 6649,
                                                                       6652, 23026, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 55692, 3, 6652,
                                                                       6655, 23032, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 55702, 3, 6655,
                                                                       6658, 23038, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 55712, 3, 6658,
                                                                       6661, 23044, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 55722, 3, 6661,
                                                                       6664, 23050, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 55732, 3, 6664,
                                                                       6667, 23056, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 55742, 3, 6667,
                                                                       6670, 23062, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 55752, 3, 6670,
                                                                       6673, 23068, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 55762, 3, 6673,
                                                                       6676, 23074, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 55772, 3, 6676,
                                                                       6679, 23080, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 55782, 3, 6679,
                                                                       6682, 23086, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 55792, 0, 3,
                                                                       55492, 22900, 55502,
                                                                       23128, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 55822, 0, 3,
                                                                       55502, 22906, 55512,
                                                                       23146, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 55852, 0, 3,
                                                                       55512, 22912, 55522,
                                                                       23164, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 55882, 0, 3,
                                                                       55522, 22918, 55532,
                                                                       23182, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 55912, 0, 3,
                                                                       55532, 22924, 55542,
                                                                       23200, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 55942, 0, 3,
                                                                       55542, 22930, 55552,
                                                                       23218, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 55972, 0, 3,
                                                                       55552, 22936, 55562,
                                                                       23236, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 56002, 0, 3,
                                                                       55562, 22942, 55572,
                                                                       23254, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 56032, 0, 3,
                                                                       55572, 22948, 55582,
                                                                       23272, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 56062, 0, 3,
                                                                       55582, 22954, 55592,
                                                                       23290, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 56092, 0, 3,
                                                                       55592, 22960, 55602,
                                                                       23308, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 56122, 0, 3,
                                                                       55602, 22966, 55612,
                                                                       23326, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 56152, 0, 3,
                                                                       55612, 22972, 55622,
                                                                       23344, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 56182, 0, 3,
                                                                       55622, 22978, 55632,
                                                                       23362, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 56212, 0, 3,
                                                                       55642, 23002, 55652,
                                                                       23416, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 56242, 0, 3,
                                                                       55652, 23008, 55662,
                                                                       23434, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 56272, 0, 3,
                                                                       55662, 23014, 55672,
                                                                       23452, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 56302, 0, 3,
                                                                       55672, 23020, 55682,
                                                                       23470, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 56332, 0, 3,
                                                                       55682, 23026, 55692,
                                                                       23488, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 56362, 0, 3,
                                                                       55692, 23032, 55702,
                                                                       23506, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 56392, 0, 3,
                                                                       55702, 23038, 55712,
                                                                       23524, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 56422, 0, 3,
                                                                       55712, 23044, 55722,
                                                                       23542, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 56452, 0, 3,
                                                                       55722, 23050, 55732,
                                                                       23560, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 56482, 0, 3,
                                                                       55732, 23056, 55742,
                                                                       23578, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 56512, 0, 3,
                                                                       55742, 23062, 55752,
                                                                       23596, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 56542, 0, 3,
                                                                       55752, 23068, 55762,
                                                                       23614, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 56572, 0, 3,
                                                                       55762, 23074, 55772,
                                                                       23632, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 56602, 0, 3,
                                                                       55772, 23080, 55782,
                                                                       23650, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 56632, 0, 3,
                                                                       55792, 23128, 55822, 6976,
                                                                       6994, 23740, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 56692, 0, 3,
                                                                       55822, 23146, 55852, 6994,
                                                                       7012, 23776, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 56752, 0, 3,
                                                                       55852, 23164, 55882, 7012,
                                                                       7030, 23812, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 56812, 0, 3,
                                                                       55882, 23182, 55912, 7030,
                                                                       7048, 23848, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 56872, 0, 3,
                                                                       55912, 23200, 55942, 7048,
                                                                       7066, 23884, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 56932, 0, 3,
                                                                       55942, 23218, 55972, 7066,
                                                                       7084, 23920, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 56992, 0, 3,
                                                                       55972, 23236, 56002, 7084,
                                                                       7102, 23956, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 57052, 0, 3,
                                                                       56002, 23254, 56032, 7102,
                                                                       7120, 23992, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 57112, 0, 3,
                                                                       56032, 23272, 56062, 7120,
                                                                       7138, 24028, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 57172, 0, 3,
                                                                       56062, 23290, 56092, 7138,
                                                                       7156, 24064, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 57232, 0, 3,
                                                                       56092, 23308, 56122, 7156,
                                                                       7174, 24100, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 57292, 0, 3,
                                                                       56122, 23326, 56152, 7174,
                                                                       7192, 24136, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 57352, 0, 3,
                                                                       56152, 23344, 56182, 7192,
                                                                       7210, 24172, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 57412, 0, 3,
                                                                       56212, 23416, 56242, 7246,
                                                                       7264, 24280, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 57472, 0, 3,
                                                                       56242, 23434, 56272, 7264,
                                                                       7282, 24316, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 57532, 0, 3,
                                                                       56272, 23452, 56302, 7282,
                                                                       7300, 24352, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 57592, 0, 3,
                                                                       56302, 23470, 56332, 7300,
                                                                       7318, 24388, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 57652, 0, 3,
                                                                       56332, 23488, 56362, 7318,
                                                                       7336, 24424, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 57712, 0, 3,
                                                                       56362, 23506, 56392, 7336,
                                                                       7354, 24460, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 57772, 0, 3,
                                                                       56392, 23524, 56422, 7354,
                                                                       7372, 24496, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 57832, 0, 3,
                                                                       56422, 23542, 56452, 7372,
                                                                       7390, 24532, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 57892, 0, 3,
                                                                       56452, 23560, 56482, 7390,
                                                                       7408, 24568, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 57952, 0, 3,
                                                                       56482, 23578, 56512, 7408,
                                                                       7426, 24604, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 58012, 0, 3,
                                                                       56512, 23596, 56542, 7426,
                                                                       7444, 24640, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 58072, 0, 3,
                                                                       56542, 23614, 56572, 7444,
                                                                       7462, 24676, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 58132, 0, 3,
                                                                       56572, 23632, 56602, 7462,
                                                                       7480, 24712, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 58192, 0, 3,
                                                                       56632, 23740, 56692, 7516,
                                                                       7546, 24868, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 58292, 0, 3,
                                                                       56692, 23776, 56752, 7546,
                                                                       7576, 24928, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 58392, 0, 3,
                                                                       56752, 23812, 56812, 7576,
                                                                       7606, 24988, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 58492, 0, 3,
                                                                       56812, 23848, 56872, 7606,
                                                                       7636, 25048, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 58592, 0, 3,
                                                                       56872, 23884, 56932, 7636,
                                                                       7666, 25108, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 58692, 0, 3,
                                                                       56932, 23920, 56992, 7666,
                                                                       7696, 25168, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 58792, 0, 3,
                                                                       56992, 23956, 57052, 7696,
                                                                       7726, 25228, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 58892, 0, 3,
                                                                       57052, 23992, 57112, 7726,
                                                                       7756, 25288, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 58992, 0, 3,
                                                                       57112, 24028, 57172, 7756,
                                                                       7786, 25348, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 59092, 0, 3,
                                                                       57172, 24064, 57232, 7786,
                                                                       7816, 25408, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 59192, 0, 3,
                                                                       57232, 24100, 57292, 7816,
                                                                       7846, 25468, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 59292, 0, 3,
                                                                       57292, 24136, 57352, 7846,
                                                                       7876, 25528, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 59392, 0, 3,
                                                                       57412, 24280, 57472, 7936,
                                                                       7966, 25708, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 59492, 0, 3,
                                                                       57472, 24316, 57532, 7966,
                                                                       7996, 25768, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 59592, 0, 3,
                                                                       57532, 24352, 57592, 7996,
                                                                       8026, 25828, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 59692, 0, 3,
                                                                       57592, 24388, 57652, 8026,
                                                                       8056, 25888, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 59792, 0, 3,
                                                                       57652, 24424, 57712, 8056,
                                                                       8086, 25948, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 59892, 0, 3,
                                                                       57712, 24460, 57772, 8086,
                                                                       8116, 26008, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 59992, 0, 3,
                                                                       57772, 24496, 57832, 8116,
                                                                       8146, 26068, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 60092, 0, 3,
                                                                       57832, 24532, 57892, 8146,
                                                                       8176, 26128, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 60192, 0, 3,
                                                                       57892, 24568, 57952, 8176,
                                                                       8206, 26188, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 60292, 0, 3,
                                                                       57952, 24604, 58012, 8206,
                                                                       8236, 26248, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 60392, 0, 3,
                                                                       58012, 24640, 58072, 8236,
                                                                       8266, 26308, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 60492, 0, 3,
                                                                       58072, 24676, 58132, 8266,
                                                                       8296, 26368, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 60592, 0, 3,
                                                                       58192, 24868, 58292, 8356,
                                                                       8401, 26608, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 60742, 0, 3,
                                                                       58292, 24928, 58392, 8401,
                                                                       8446, 26698, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 60892, 0, 3,
                                                                       58392, 24988, 58492, 8446,
                                                                       8491, 26788, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 61042, 0, 3,
                                                                       58492, 25048, 58592, 8491,
                                                                       8536, 26878, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 61192, 0, 3,
                                                                       58592, 25108, 58692, 8536,
                                                                       8581, 26968, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 61342, 0, 3,
                                                                       58692, 25168, 58792, 8581,
                                                                       8626, 27058, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 61492, 0, 3,
                                                                       58792, 25228, 58892, 8626,
                                                                       8671, 27148, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 61642, 0, 3,
                                                                       58892, 25288, 58992, 8671,
                                                                       8716, 27238, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 61792, 0, 3,
                                                                       58992, 25348, 59092, 8716,
                                                                       8761, 27328, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 61942, 0, 3,
                                                                       59092, 25408, 59192, 8761,
                                                                       8806, 27418, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 62092, 0, 3,
                                                                       59192, 25468, 59292, 8806,
                                                                       8851, 27508, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 62242, 0, 3,
                                                                       59392, 25708, 59492, 8941,
                                                                       8986, 27778, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 62392, 0, 3,
                                                                       59492, 25768, 59592, 8986,
                                                                       9031, 27868, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 62542, 0, 3,
                                                                       59592, 25828, 59692, 9031,
                                                                       9076, 27958, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 62692, 0, 3,
                                                                       59692, 25888, 59792, 9076,
                                                                       9121, 28048, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 62842, 0, 3,
                                                                       59792, 25948, 59892, 9121,
                                                                       9166, 28138, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 62992, 0, 3,
                                                                       59892, 26008, 59992, 9166,
                                                                       9211, 28228, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 63142, 0, 3,
                                                                       59992, 26068, 60092, 9211,
                                                                       9256, 28318, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 63292, 0, 3,
                                                                       60092, 26128, 60192, 9256,
                                                                       9301, 28408, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 63442, 0, 3,
                                                                       60192, 26188, 60292, 9301,
                                                                       9346, 28498, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 63592, 0, 3,
                                                                       60292, 26248, 60392, 9346,
                                                                       9391, 28588, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 63742, 0, 3,
                                                                       60392, 26308, 60492, 9391,
                                                                       9436, 28678, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 63892, 0, 3,
                                                                       60592, 26608, 60742, 9526,
                                                                       9589, 29020, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 64102, 0, 3,
                                                                       60742, 26698, 60892, 9589,
                                                                       9652, 29146, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 64312, 0, 3,
                                                                       60892, 26788, 61042, 9652,
                                                                       9715, 29272, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 64522, 0, 3,
                                                                       61042, 26878, 61192, 9715,
                                                                       9778, 29398, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 64732, 0, 3,
                                                                       61192, 26968, 61342, 9778,
                                                                       9841, 29524, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 64942, 0, 3,
                                                                       61342, 27058, 61492, 9841,
                                                                       9904, 29650, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 65152, 0, 3,
                                                                       61492, 27148, 61642, 9904,
                                                                       9967, 29776, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 65362, 0, 3,
                                                                       61642, 27238, 61792, 9967,
                                                                       10030, 29902, ncols,
                                                                       gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 65572, 0, 3,
                                                                       61792, 27328, 61942,
                                                                       10030, 10093, 30028,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 65782, 0, 3,
                                                                       61942, 27418, 62092,
                                                                       10093, 10156, 30154,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 65992, 0, 3,
                                                                       62242, 27778, 62392,
                                                                       10282, 10345, 30532,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 66202, 0, 3,
                                                                       62392, 27868, 62542,
                                                                       10345, 10408, 30658,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 66412, 0, 3,
                                                                       62542, 27958, 62692,
                                                                       10408, 10471, 30784,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 66622, 0, 3,
                                                                       62692, 28048, 62842,
                                                                       10471, 10534, 30910,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 66832, 0, 3,
                                                                       62842, 28138, 62992,
                                                                       10534, 10597, 31036,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 67042, 0, 3,
                                                                       62992, 28228, 63142,
                                                                       10597, 10660, 31162,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 67252, 0, 3,
                                                                       63142, 28318, 63292,
                                                                       10660, 10723, 31288,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 67462, 0, 3,
                                                                       63292, 28408, 63442,
                                                                       10723, 10786, 31414,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 67672, 0, 3,
                                                                       63442, 28498, 63592,
                                                                       10786, 10849, 31540,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 67882, 0, 3,
                                                                       63592, 28588, 63742,
                                                                       10849, 10912, 31666,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 68092, 0, 3,
                                                                       63892, 29020, 64102,
                                                                       11038, 11122, 32128,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 68372, 0, 3,
                                                                       64102, 29146, 64312,
                                                                       11122, 11206, 32296,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 68652, 0, 3,
                                                                       64312, 29272, 64522,
                                                                       11206, 11290, 32464,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 68932, 0, 3,
                                                                       64522, 29398, 64732,
                                                                       11290, 11374, 32632,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 69212, 0, 3,
                                                                       64732, 29524, 64942,
                                                                       11374, 11458, 32800,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 69492, 0, 3,
                                                                       64942, 29650, 65152,
                                                                       11458, 11542, 32968,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 69772, 0, 3,
                                                                       65152, 29776, 65362,
                                                                       11542, 11626, 33136,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 70052, 0, 3,
                                                                       65362, 29902, 65572,
                                                                       11626, 11710, 33304,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 70332, 0, 3,
                                                                       65572, 30028, 65782,
                                                                       11710, 11794, 33472,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 70612, 0, 3,
                                                                       65992, 30532, 66202,
                                                                       11962, 12046, 33976,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 70892, 0, 3,
                                                                       66202, 30658, 66412,
                                                                       12046, 12130, 34144,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 71172, 0, 3,
                                                                       66412, 30784, 66622,
                                                                       12130, 12214, 34312,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 71452, 0, 3,
                                                                       66622, 30910, 66832,
                                                                       12214, 12298, 34480,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 71732, 0, 3,
                                                                       66832, 31036, 67042,
                                                                       12298, 12382, 34648,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 72012, 0, 3,
                                                                       67042, 31162, 67252,
                                                                       12382, 12466, 34816,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 72292, 0, 3,
                                                                       67252, 31288, 67462,
                                                                       12466, 12550, 34984,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 72572, 0, 3,
                                                                       67462, 31414, 67672,
                                                                       12550, 12634, 35152,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 72852, 0, 3,
                                                                       67672, 31540, 67882,
                                                                       12634, 12718, 35320,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 73132, 0, 3,
                                                                       68092, 32128, 68372,
                                                                       12886, 12994, 35920,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 73492, 0, 3,
                                                                       68372, 32296, 68652,
                                                                       12994, 13102, 36136,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 73852, 0, 3,
                                                                       68652, 32464, 68932,
                                                                       13102, 13210, 36352,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 74212, 0, 3,
                                                                       68932, 32632, 69212,
                                                                       13210, 13318, 36568,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 74572, 0, 3,
                                                                       69212, 32800, 69492,
                                                                       13318, 13426, 36784,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 74932, 0, 3,
                                                                       69492, 32968, 69772,
                                                                       13426, 13534, 37000,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 75292, 0, 3,
                                                                       69772, 33136, 70052,
                                                                       13534, 13642, 37216,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 75652, 0, 3,
                                                                       70052, 33304, 70332,
                                                                       13642, 13750, 37432,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 76012, 0, 3,
                                                                       70612, 33976, 70892,
                                                                       13966, 14074, 38080,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 76372, 0, 3,
                                                                       70892, 34144, 71172,
                                                                       14074, 14182, 38296,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 76732, 0, 3,
                                                                       71172, 34312, 71452,
                                                                       14182, 14290, 38512,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 77092, 0, 3,
                                                                       71452, 34480, 71732,
                                                                       14290, 14398, 38728,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 77452, 0, 3,
                                                                       71732, 34648, 72012,
                                                                       14398, 14506, 38944,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 77812, 0, 3,
                                                                       72012, 34816, 72292,
                                                                       14506, 14614, 39160,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 78172, 0, 3,
                                                                       72292, 34984, 72572,
                                                                       14614, 14722, 39376,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 78532, 0, 3,
                                                                       72572, 35152, 72852,
                                                                       14722, 14830, 39592,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 78892, 0, 3,
                                                                       73132, 35920, 73492,
                                                                       15046, 15181, 40348,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 79342, 0, 3,
                                                                       73492, 36136, 73852,
                                                                       15181, 15316, 40618,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 79792, 0, 3,
                                                                       73852, 36352, 74212,
                                                                       15316, 15451, 40888,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 80242, 0, 3,
                                                                       74212, 36568, 74572,
                                                                       15451, 15586, 41158,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 80692, 0, 3,
                                                                       74572, 36784, 74932,
                                                                       15586, 15721, 41428,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 81142, 0, 3,
                                                                       74932, 37000, 75292,
                                                                       15721, 15856, 41698,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 81592, 0, 3,
                                                                       75292, 37216, 75652,
                                                                       15856, 15991, 41968,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 82042, 0, 3,
                                                                       76012, 38080, 76372,
                                                                       16261, 16396, 42778,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 82492, 0, 3,
                                                                       76372, 38296, 76732,
                                                                       16396, 16531, 43048,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 82942, 0, 3,
                                                                       76732, 38512, 77092,
                                                                       16531, 16666, 43318,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 83392, 0, 3,
                                                                       77092, 38728, 77452,
                                                                       16666, 16801, 43588,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 83842, 0, 3,
                                                                       77452, 38944, 77812,
                                                                       16801, 16936, 43858,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 84292, 0, 3,
                                                                       77812, 39160, 78172,
                                                                       16936, 17071, 44128,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 84742, 0, 3,
                                                                       78172, 39376, 78532,
                                                                       17071, 17206, 44398,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 85192, 0, 3,
                                                                       78892, 40348, 79342,
                                                                       17476, 17641, 45328,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 85742, 0, 3,
                                                                       79342, 40618, 79792,
                                                                       17641, 17806, 45658,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 86292, 0, 3,
                                                                       79792, 40888, 80242,
                                                                       17806, 17971, 45988,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 86842, 0, 3,
                                                                       80242, 41158, 80692,
                                                                       17971, 18136, 46318,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 87392, 0, 3,
                                                                       80692, 41428, 81142,
                                                                       18136, 18301, 46648,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 87942, 0, 3,
                                                                       81142, 41698, 81592,
                                                                       18301, 18466, 46978,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 88492, 0, 3,
                                                                       82042, 42778, 82492,
                                                                       18796, 18961, 47968,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 89042, 0, 3,
                                                                       82492, 43048, 82942,
                                                                       18961, 19126, 48298,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 89592, 0, 3,
                                                                       82942, 43318, 83392,
                                                                       19126, 19291, 48628,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 90142, 0, 3,
                                                                       83392, 43588, 83842,
                                                                       19291, 19456, 48958,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 90692, 0, 3,
                                                                       83842, 43858, 84292,
                                                                       19456, 19621, 49288,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 91242, 0, 3,
                                                                       84292, 44128, 84742,
                                                                       19621, 19786, 49618,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 91792, 0, 3,
                                                                       85192, 45328, 85742,
                                                                       20116, 20314, 50740,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 92452, 0, 3,
                                                                       85742, 45658, 86292,
                                                                       20314, 20512, 51136,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 93112, 0, 3,
                                                                       86292, 45988, 86842,
                                                                       20512, 20710, 51532,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 93772, 0, 3,
                                                                       86842, 46318, 87392,
                                                                       20710, 20908, 51928,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 94432, 0, 3,
                                                                       87392, 46648, 87942,
                                                                       20908, 21106, 52324,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 95092, 0, 3,
                                                                       88492, 47968, 89042,
                                                                       21502, 21700, 53512,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 95752, 0, 3,
                                                                       89042, 48298, 89592,
                                                                       21700, 21898, 53908,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 96412, 0, 3,
                                                                       89592, 48628, 90142,
                                                                       21898, 22096, 54304,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 97072, 0, 3,
                                                                       90142, 48958, 90692,
                                                                       22096, 22294, 54700,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 97732, 0, 3,
                                                                       90692, 49288, 91242,
                                                                       22294, 22492, 55096,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 98392, 3, 22888,
                                                                       22894, 55492, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 98407, 3, 22894,
                                                                       22900, 55502, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 98422, 3, 22900,
                                                                       22906, 55512, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 98437, 3, 22906,
                                                                       22912, 55522, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 98452, 3, 22912,
                                                                       22918, 55532, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 98467, 3, 22918,
                                                                       22924, 55542, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 98482, 3, 22924,
                                                                       22930, 55552, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 98497, 3, 22930,
                                                                       22936, 55562, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 98512, 3, 22936,
                                                                       22942, 55572, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 98527, 3, 22942,
                                                                       22948, 55582, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 98542, 3, 22948,
                                                                       22954, 55592, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 98557, 3, 22954,
                                                                       22960, 55602, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 98572, 3, 22960,
                                                                       22966, 55612, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 98587, 3, 22966,
                                                                       22972, 55622, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 98602, 3, 22972,
                                                                       22978, 55632, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 98617, 3, 22990,
                                                                       22996, 55642, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 98632, 3, 22996,
                                                                       23002, 55652, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 98647, 3, 23002,
                                                                       23008, 55662, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 98662, 3, 23008,
                                                                       23014, 55672, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 98677, 3, 23014,
                                                                       23020, 55682, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 98692, 3, 23020,
                                                                       23026, 55692, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 98707, 3, 23026,
                                                                       23032, 55702, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 98722, 3, 23032,
                                                                       23038, 55712, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 98737, 3, 23038,
                                                                       23044, 55722, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 98752, 3, 23044,
                                                                       23050, 55732, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 98767, 3, 23050,
                                                                       23056, 55742, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 98782, 3, 23056,
                                                                       23062, 55752, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 98797, 3, 23062,
                                                                       23068, 55762, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 98812, 3, 23068,
                                                                       23074, 55772, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 98827, 3, 23074,
                                                                       23080, 55782, ncols,
                                                                       gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 98842, 0, 3,
                                                                       98392, 55492, 98407,
                                                                       23092, 23110, 55792,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 98887, 0, 3,
                                                                       98407, 55502, 98422,
                                                                       23110, 23128, 55822,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 98932, 0, 3,
                                                                       98422, 55512, 98437,
                                                                       23128, 23146, 55852,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 98977, 0, 3,
                                                                       98437, 55522, 98452,
                                                                       23146, 23164, 55882,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 99022, 0, 3,
                                                                       98452, 55532, 98467,
                                                                       23164, 23182, 55912,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 99067, 0, 3,
                                                                       98467, 55542, 98482,
                                                                       23182, 23200, 55942,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 99112, 0, 3,
                                                                       98482, 55552, 98497,
                                                                       23200, 23218, 55972,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 99157, 0, 3,
                                                                       98497, 55562, 98512,
                                                                       23218, 23236, 56002,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 99202, 0, 3,
                                                                       98512, 55572, 98527,
                                                                       23236, 23254, 56032,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 99247, 0, 3,
                                                                       98527, 55582, 98542,
                                                                       23254, 23272, 56062,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 99292, 0, 3,
                                                                       98542, 55592, 98557,
                                                                       23272, 23290, 56092,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 99337, 0, 3,
                                                                       98557, 55602, 98572,
                                                                       23290, 23308, 56122,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 99382, 0, 3,
                                                                       98572, 55612, 98587,
                                                                       23308, 23326, 56152,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 99427, 0, 3,
                                                                       98587, 55622, 98602,
                                                                       23326, 23344, 56182,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 99472, 0, 3,
                                                                       98617, 55642, 98632,
                                                                       23380, 23398, 56212,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 99517, 0, 3,
                                                                       98632, 55652, 98647,
                                                                       23398, 23416, 56242,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 99562, 0, 3,
                                                                       98647, 55662, 98662,
                                                                       23416, 23434, 56272,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 99607, 0, 3,
                                                                       98662, 55672, 98677,
                                                                       23434, 23452, 56302,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 99652, 0, 3,
                                                                       98677, 55682, 98692,
                                                                       23452, 23470, 56332,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 99697, 0, 3,
                                                                       98692, 55692, 98707,
                                                                       23470, 23488, 56362,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 99742, 0, 3,
                                                                       98707, 55702, 98722,
                                                                       23488, 23506, 56392,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 99787, 0, 3,
                                                                       98722, 55712, 98737,
                                                                       23506, 23524, 56422,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 99832, 0, 3,
                                                                       98737, 55722, 98752,
                                                                       23524, 23542, 56452,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 99877, 0, 3,
                                                                       98752, 55732, 98767,
                                                                       23542, 23560, 56482,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 99922, 0, 3,
                                                                       98767, 55742, 98782,
                                                                       23560, 23578, 56512,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 99967, 0, 3,
                                                                       98782, 55752, 98797,
                                                                       23578, 23596, 56542,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 100012, 0, 3,
                                                                       98797, 55762, 98812,
                                                                       23596, 23614, 56572,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 100057, 0, 3,
                                                                       98812, 55772, 98827,
                                                                       23614, 23632, 56602,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 100102, 0, 3,
                                                                       98842, 55792, 98887,
                                                                       23668, 23704, 56632,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 100192, 0, 3,
                                                                       98887, 55822, 98932,
                                                                       23704, 23740, 56692,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 100282, 0, 3,
                                                                       98932, 55852, 98977,
                                                                       23740, 23776, 56752,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 100372, 0, 3,
                                                                       98977, 55882, 99022,
                                                                       23776, 23812, 56812,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 100462, 0, 3,
                                                                       99022, 55912, 99067,
                                                                       23812, 23848, 56872,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 100552, 0, 3,
                                                                       99067, 55942, 99112,
                                                                       23848, 23884, 56932,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 100642, 0, 3,
                                                                       99112, 55972, 99157,
                                                                       23884, 23920, 56992,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 100732, 0, 3,
                                                                       99157, 56002, 99202,
                                                                       23920, 23956, 57052,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 100822, 0, 3,
                                                                       99202, 56032, 99247,
                                                                       23956, 23992, 57112,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 100912, 0, 3,
                                                                       99247, 56062, 99292,
                                                                       23992, 24028, 57172,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 101002, 0, 3,
                                                                       99292, 56092, 99337,
                                                                       24028, 24064, 57232,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 101092, 0, 3,
                                                                       99337, 56122, 99382,
                                                                       24064, 24100, 57292,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 101182, 0, 3,
                                                                       99382, 56152, 99427,
                                                                       24100, 24136, 57352,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 101272, 0, 3,
                                                                       99472, 56212, 99517,
                                                                       24208, 24244, 57412,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 101362, 0, 3,
                                                                       99517, 56242, 99562,
                                                                       24244, 24280, 57472,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 101452, 0, 3,
                                                                       99562, 56272, 99607,
                                                                       24280, 24316, 57532,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 101542, 0, 3,
                                                                       99607, 56302, 99652,
                                                                       24316, 24352, 57592,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 101632, 0, 3,
                                                                       99652, 56332, 99697,
                                                                       24352, 24388, 57652,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 101722, 0, 3,
                                                                       99697, 56362, 99742,
                                                                       24388, 24424, 57712,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 101812, 0, 3,
                                                                       99742, 56392, 99787,
                                                                       24424, 24460, 57772,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 101902, 0, 3,
                                                                       99787, 56422, 99832,
                                                                       24460, 24496, 57832,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 101992, 0, 3,
                                                                       99832, 56452, 99877,
                                                                       24496, 24532, 57892,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 102082, 0, 3,
                                                                       99877, 56482, 99922,
                                                                       24532, 24568, 57952,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 102172, 0, 3,
                                                                       99922, 56512, 99967,
                                                                       24568, 24604, 58012,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 102262, 0, 3,
                                                                       99967, 56542, 100012,
                                                                       24604, 24640, 58072,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 102352, 0, 3,
                                                                       100012, 56572, 100057,
                                                                       24640, 24676, 58132,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 102442, 0, 3,
                                                                       100102, 56632, 100192,
                                                                       24748, 24808, 58192,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 102592, 0, 3,
                                                                       100192, 56692, 100282,
                                                                       24808, 24868, 58292,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 102742, 0, 3,
                                                                       100282, 56752, 100372,
                                                                       24868, 24928, 58392,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 102892, 0, 3,
                                                                       100372, 56812, 100462,
                                                                       24928, 24988, 58492,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 103042, 0, 3,
                                                                       100462, 56872, 100552,
                                                                       24988, 25048, 58592,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 103192, 0, 3,
                                                                       100552, 56932, 100642,
                                                                       25048, 25108, 58692,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 103342, 0, 3,
                                                                       100642, 56992, 100732,
                                                                       25108, 25168, 58792,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 103492, 0, 3,
                                                                       100732, 57052, 100822,
                                                                       25168, 25228, 58892,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 103642, 0, 3,
                                                                       100822, 57112, 100912,
                                                                       25228, 25288, 58992,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 103792, 0, 3,
                                                                       100912, 57172, 101002,
                                                                       25288, 25348, 59092,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 103942, 0, 3,
                                                                       101002, 57232, 101092,
                                                                       25348, 25408, 59192,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 104092, 0, 3,
                                                                       101092, 57292, 101182,
                                                                       25408, 25468, 59292,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 104242, 0, 3,
                                                                       101272, 57412, 101362,
                                                                       25588, 25648, 59392,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 104392, 0, 3,
                                                                       101362, 57472, 101452,
                                                                       25648, 25708, 59492,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 104542, 0, 3,
                                                                       101452, 57532, 101542,
                                                                       25708, 25768, 59592,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 104692, 0, 3,
                                                                       101542, 57592, 101632,
                                                                       25768, 25828, 59692,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 104842, 0, 3,
                                                                       101632, 57652, 101722,
                                                                       25828, 25888, 59792,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 104992, 0, 3,
                                                                       101722, 57712, 101812,
                                                                       25888, 25948, 59892,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 105142, 0, 3,
                                                                       101812, 57772, 101902,
                                                                       25948, 26008, 59992,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 105292, 0, 3,
                                                                       101902, 57832, 101992,
                                                                       26008, 26068, 60092,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 105442, 0, 3,
                                                                       101992, 57892, 102082,
                                                                       26068, 26128, 60192,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 105592, 0, 3,
                                                                       102082, 57952, 102172,
                                                                       26128, 26188, 60292,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 105742, 0, 3,
                                                                       102172, 58012, 102262,
                                                                       26188, 26248, 60392,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 105892, 0, 3,
                                                                       102262, 58072, 102352,
                                                                       26248, 26308, 60492,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 106042, 0, 3,
                                                                       102442, 58192, 102592,
                                                                       26428, 26518, 60592,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 106267, 0, 3,
                                                                       102592, 58292, 102742,
                                                                       26518, 26608, 60742,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 106492, 0, 3,
                                                                       102742, 58392, 102892,
                                                                       26608, 26698, 60892,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 106717, 0, 3,
                                                                       102892, 58492, 103042,
                                                                       26698, 26788, 61042,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 106942, 0, 3,
                                                                       103042, 58592, 103192,
                                                                       26788, 26878, 61192,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 107167, 0, 3,
                                                                       103192, 58692, 103342,
                                                                       26878, 26968, 61342,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 107392, 0, 3,
                                                                       103342, 58792, 103492,
                                                                       26968, 27058, 61492,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 107617, 0, 3,
                                                                       103492, 58892, 103642,
                                                                       27058, 27148, 61642,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 107842, 0, 3,
                                                                       103642, 58992, 103792,
                                                                       27148, 27238, 61792,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 108067, 0, 3,
                                                                       103792, 59092, 103942,
                                                                       27238, 27328, 61942,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 108292, 0, 3,
                                                                       103942, 59192, 104092,
                                                                       27328, 27418, 62092,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 108517, 0, 3,
                                                                       104242, 59392, 104392,
                                                                       27598, 27688, 62242,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 108742, 0, 3,
                                                                       104392, 59492, 104542,
                                                                       27688, 27778, 62392,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 108967, 0, 3,
                                                                       104542, 59592, 104692,
                                                                       27778, 27868, 62542,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 109192, 0, 3,
                                                                       104692, 59692, 104842,
                                                                       27868, 27958, 62692,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 109417, 0, 3,
                                                                       104842, 59792, 104992,
                                                                       27958, 28048, 62842,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 109642, 0, 3,
                                                                       104992, 59892, 105142,
                                                                       28048, 28138, 62992,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 109867, 0, 3,
                                                                       105142, 59992, 105292,
                                                                       28138, 28228, 63142,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 110092, 0, 3,
                                                                       105292, 60092, 105442,
                                                                       28228, 28318, 63292,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 110317, 0, 3,
                                                                       105442, 60192, 105592,
                                                                       28318, 28408, 63442,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 110542, 0, 3,
                                                                       105592, 60292, 105742,
                                                                       28408, 28498, 63592,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 110767, 0, 3,
                                                                       105742, 60392, 105892,
                                                                       28498, 28588, 63742,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 110992, 0, 3,
                                                                       106042, 60592, 106267,
                                                                       28768, 28894, 63892,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 111307, 0, 3,
                                                                       106267, 60742, 106492,
                                                                       28894, 29020, 64102,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 111622, 0, 3,
                                                                       106492, 60892, 106717,
                                                                       29020, 29146, 64312,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 111937, 0, 3,
                                                                       106717, 61042, 106942,
                                                                       29146, 29272, 64522,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 112252, 0, 3,
                                                                       106942, 61192, 107167,
                                                                       29272, 29398, 64732,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 112567, 0, 3,
                                                                       107167, 61342, 107392,
                                                                       29398, 29524, 64942,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 112882, 0, 3,
                                                                       107392, 61492, 107617,
                                                                       29524, 29650, 65152,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 113197, 0, 3,
                                                                       107617, 61642, 107842,
                                                                       29650, 29776, 65362,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 113512, 0, 3,
                                                                       107842, 61792, 108067,
                                                                       29776, 29902, 65572,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 113827, 0, 3,
                                                                       108067, 61942, 108292,
                                                                       29902, 30028, 65782,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 114142, 0, 3,
                                                                       108517, 62242, 108742,
                                                                       30280, 30406, 65992,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 114457, 0, 3,
                                                                       108742, 62392, 108967,
                                                                       30406, 30532, 66202,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 114772, 0, 3,
                                                                       108967, 62542, 109192,
                                                                       30532, 30658, 66412,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 115087, 0, 3,
                                                                       109192, 62692, 109417,
                                                                       30658, 30784, 66622,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 115402, 0, 3,
                                                                       109417, 62842, 109642,
                                                                       30784, 30910, 66832,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 115717, 0, 3,
                                                                       109642, 62992, 109867,
                                                                       30910, 31036, 67042,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 116032, 0, 3,
                                                                       109867, 63142, 110092,
                                                                       31036, 31162, 67252,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 116347, 0, 3,
                                                                       110092, 63292, 110317,
                                                                       31162, 31288, 67462,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 116662, 0, 3,
                                                                       110317, 63442, 110542,
                                                                       31288, 31414, 67672,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 116977, 0, 3,
                                                                       110542, 63592, 110767,
                                                                       31414, 31540, 67882,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 117292, 0, 3,
                                                                       110992, 63892, 111307,
                                                                       31792, 31960, 68092,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 117712, 0, 3,
                                                                       111307, 64102, 111622,
                                                                       31960, 32128, 68372,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 118132, 0, 3,
                                                                       111622, 64312, 111937,
                                                                       32128, 32296, 68652,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 118552, 0, 3,
                                                                       111937, 64522, 112252,
                                                                       32296, 32464, 68932,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 118972, 0, 3,
                                                                       112252, 64732, 112567,
                                                                       32464, 32632, 69212,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 119392, 0, 3,
                                                                       112567, 64942, 112882,
                                                                       32632, 32800, 69492,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 119812, 0, 3,
                                                                       112882, 65152, 113197,
                                                                       32800, 32968, 69772,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 120232, 0, 3,
                                                                       113197, 65362, 113512,
                                                                       32968, 33136, 70052,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 120652, 0, 3,
                                                                       113512, 65572, 113827,
                                                                       33136, 33304, 70332,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 121072, 0, 3,
                                                                       114142, 65992, 114457,
                                                                       33640, 33808, 70612,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 121492, 0, 3,
                                                                       114457, 66202, 114772,
                                                                       33808, 33976, 70892,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 121912, 0, 3,
                                                                       114772, 66412, 115087,
                                                                       33976, 34144, 71172,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 122332, 0, 3,
                                                                       115087, 66622, 115402,
                                                                       34144, 34312, 71452,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 122752, 0, 3,
                                                                       115402, 66832, 115717,
                                                                       34312, 34480, 71732,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 123172, 0, 3,
                                                                       115717, 67042, 116032,
                                                                       34480, 34648, 72012,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 123592, 0, 3,
                                                                       116032, 67252, 116347,
                                                                       34648, 34816, 72292,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 124012, 0, 3,
                                                                       116347, 67462, 116662,
                                                                       34816, 34984, 72572,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 124432, 0, 3,
                                                                       116662, 67672, 116977,
                                                                       34984, 35152, 72852,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 124852, 0, 3,
                                                                       117292, 68092, 117712,
                                                                       35488, 35704, 73132,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 125392, 0, 3,
                                                                       117712, 68372, 118132,
                                                                       35704, 35920, 73492,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 125932, 0, 3,
                                                                       118132, 68652, 118552,
                                                                       35920, 36136, 73852,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 126472, 0, 3,
                                                                       118552, 68932, 118972,
                                                                       36136, 36352, 74212,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 127012, 0, 3,
                                                                       118972, 69212, 119392,
                                                                       36352, 36568, 74572,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 127552, 0, 3,
                                                                       119392, 69492, 119812,
                                                                       36568, 36784, 74932,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 128092, 0, 3,
                                                                       119812, 69772, 120232,
                                                                       36784, 37000, 75292,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 128632, 0, 3,
                                                                       120232, 70052, 120652,
                                                                       37000, 37216, 75652,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 129172, 0, 3,
                                                                       121072, 70612, 121492,
                                                                       37648, 37864, 76012,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 129712, 0, 3,
                                                                       121492, 70892, 121912,
                                                                       37864, 38080, 76372,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 130252, 0, 3,
                                                                       121912, 71172, 122332,
                                                                       38080, 38296, 76732,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 130792, 0, 3,
                                                                       122332, 71452, 122752,
                                                                       38296, 38512, 77092,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 131332, 0, 3,
                                                                       122752, 71732, 123172,
                                                                       38512, 38728, 77452,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 131872, 0, 3,
                                                                       123172, 72012, 123592,
                                                                       38728, 38944, 77812,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 132412, 0, 3,
                                                                       123592, 72292, 124012,
                                                                       38944, 39160, 78172,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 132952, 0, 3,
                                                                       124012, 72572, 124432,
                                                                       39160, 39376, 78532,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 133492, 0, 3,
                                                                       124852, 73132, 125392,
                                                                       39808, 40078, 78892,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 134167, 0, 3,
                                                                       125392, 73492, 125932,
                                                                       40078, 40348, 79342,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 134842, 0, 3,
                                                                       125932, 73852, 126472,
                                                                       40348, 40618, 79792,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 135517, 0, 3,
                                                                       126472, 74212, 127012,
                                                                       40618, 40888, 80242,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 136192, 0, 3,
                                                                       127012, 74572, 127552,
                                                                       40888, 41158, 80692,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 136867, 0, 3,
                                                                       127552, 74932, 128092,
                                                                       41158, 41428, 81142,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 137542, 0, 3,
                                                                       128092, 75292, 128632,
                                                                       41428, 41698, 81592,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 138217, 0, 3,
                                                                       129172, 76012, 129712,
                                                                       42238, 42508, 82042,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 138892, 0, 3,
                                                                       129712, 76372, 130252,
                                                                       42508, 42778, 82492,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 139567, 0, 3,
                                                                       130252, 76732, 130792,
                                                                       42778, 43048, 82942,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 140242, 0, 3,
                                                                       130792, 77092, 131332,
                                                                       43048, 43318, 83392,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 140917, 0, 3,
                                                                       131332, 77452, 131872,
                                                                       43318, 43588, 83842,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 141592, 0, 3,
                                                                       131872, 77812, 132412,
                                                                       43588, 43858, 84292,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 142267, 0, 3,
                                                                       132412, 78172, 132952,
                                                                       43858, 44128, 84742,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 142942, 0, 3,
                                                                       133492, 78892, 134167,
                                                                       44668, 44998, 85192,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 143767, 0, 3,
                                                                       134167, 79342, 134842,
                                                                       44998, 45328, 85742,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 144592, 0, 3,
                                                                       134842, 79792, 135517,
                                                                       45328, 45658, 86292,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 145417, 0, 3,
                                                                       135517, 80242, 136192,
                                                                       45658, 45988, 86842,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 146242, 0, 3,
                                                                       136192, 80692, 136867,
                                                                       45988, 46318, 87392,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 147067, 0, 3,
                                                                       136867, 81142, 137542,
                                                                       46318, 46648, 87942,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 147892, 0, 3,
                                                                       138217, 82042, 138892,
                                                                       47308, 47638, 88492,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 148717, 0, 3,
                                                                       138892, 82492, 139567,
                                                                       47638, 47968, 89042,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 149542, 0, 3,
                                                                       139567, 82942, 140242,
                                                                       47968, 48298, 89592,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 150367, 0, 3,
                                                                       140242, 83392, 140917,
                                                                       48298, 48628, 90142,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 151192, 0, 3,
                                                                       140917, 83842, 141592,
                                                                       48628, 48958, 90692,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 152017, 0, 3,
                                                                       141592, 84292, 142267,
                                                                       48958, 49288, 91242,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 152842, 0, 3,
                                                                       142942, 85192, 143767,
                                                                       49948, 50344, 91792,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 153832, 0, 3,
                                                                       143767, 85742, 144592,
                                                                       50344, 50740, 92452,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 154822, 0, 3,
                                                                       144592, 86292, 145417,
                                                                       50740, 51136, 93112,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 155812, 0, 3,
                                                                       145417, 86842, 146242,
                                                                       51136, 51532, 93772,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 156802, 0, 3,
                                                                       146242, 87392, 147067,
                                                                       51532, 51928, 94432,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 157792, 0, 3,
                                                                       147892, 88492, 148717,
                                                                       52720, 53116, 95092,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 158782, 0, 3,
                                                                       148717, 89042, 149542,
                                                                       53116, 53512, 95752,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 159772, 0, 3,
                                                                       149542, 89592, 150367,
                                                                       53512, 53908, 96412,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 160762, 0, 3,
                                                                       150367, 90142, 151192,
                                                                       53908, 54304, 97072,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 161752, 0, 3,
                                                                       151192, 90692, 152017,
                                                                       54304, 54700, 97732,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 162742, 3, 55492,
                                                                       55502, 98422, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 162763, 3, 55502,
                                                                       55512, 98437, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 162784, 3, 55512,
                                                                       55522, 98452, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 162805, 3, 55522,
                                                                       55532, 98467, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 162826, 3, 55532,
                                                                       55542, 98482, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 162847, 3, 55542,
                                                                       55552, 98497, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 162868, 3, 55552,
                                                                       55562, 98512, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 162889, 3, 55562,
                                                                       55572, 98527, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 162910, 3, 55572,
                                                                       55582, 98542, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 162931, 3, 55582,
                                                                       55592, 98557, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 162952, 3, 55592,
                                                                       55602, 98572, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 162973, 3, 55602,
                                                                       55612, 98587, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 162994, 3, 55612,
                                                                       55622, 98602, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 163015, 3, 55642,
                                                                       55652, 98647, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 163036, 3, 55652,
                                                                       55662, 98662, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 163057, 3, 55662,
                                                                       55672, 98677, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 163078, 3, 55672,
                                                                       55682, 98692, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 163099, 3, 55682,
                                                                       55692, 98707, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 163120, 3, 55692,
                                                                       55702, 98722, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 163141, 3, 55702,
                                                                       55712, 98737, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 163162, 3, 55712,
                                                                       55722, 98752, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 163183, 3, 55722,
                                                                       55732, 98767, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 163204, 3, 55732,
                                                                       55742, 98782, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 163225, 3, 55742,
                                                                       55752, 98797, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 163246, 3, 55752,
                                                                       55762, 98812, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 163267, 3, 55762,
                                                                       55772, 98827, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 163288, 0, 3,
                                                                       162742, 98422, 162763,
                                                                       55792, 55822, 98932,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 163351, 0, 3,
                                                                       162763, 98437, 162784,
                                                                       55822, 55852, 98977,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 163414, 0, 3,
                                                                       162784, 98452, 162805,
                                                                       55852, 55882, 99022,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 163477, 0, 3,
                                                                       162805, 98467, 162826,
                                                                       55882, 55912, 99067,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 163540, 0, 3,
                                                                       162826, 98482, 162847,
                                                                       55912, 55942, 99112,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 163603, 0, 3,
                                                                       162847, 98497, 162868,
                                                                       55942, 55972, 99157,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 163666, 0, 3,
                                                                       162868, 98512, 162889,
                                                                       55972, 56002, 99202,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 163729, 0, 3,
                                                                       162889, 98527, 162910,
                                                                       56002, 56032, 99247,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 163792, 0, 3,
                                                                       162910, 98542, 162931,
                                                                       56032, 56062, 99292,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 163855, 0, 3,
                                                                       162931, 98557, 162952,
                                                                       56062, 56092, 99337,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 163918, 0, 3,
                                                                       162952, 98572, 162973,
                                                                       56092, 56122, 99382,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 163981, 0, 3,
                                                                       162973, 98587, 162994,
                                                                       56122, 56152, 99427,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 164044, 0, 3,
                                                                       163015, 98647, 163036,
                                                                       56212, 56242, 99562,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 164107, 0, 3,
                                                                       163036, 98662, 163057,
                                                                       56242, 56272, 99607,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 164170, 0, 3,
                                                                       163057, 98677, 163078,
                                                                       56272, 56302, 99652,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 164233, 0, 3,
                                                                       163078, 98692, 163099,
                                                                       56302, 56332, 99697,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 164296, 0, 3,
                                                                       163099, 98707, 163120,
                                                                       56332, 56362, 99742,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 164359, 0, 3,
                                                                       163120, 98722, 163141,
                                                                       56362, 56392, 99787,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 164422, 0, 3,
                                                                       163141, 98737, 163162,
                                                                       56392, 56422, 99832,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 164485, 0, 3,
                                                                       163162, 98752, 163183,
                                                                       56422, 56452, 99877,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 164548, 0, 3,
                                                                       163183, 98767, 163204,
                                                                       56452, 56482, 99922,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 164611, 0, 3,
                                                                       163204, 98782, 163225,
                                                                       56482, 56512, 99967,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 164674, 0, 3,
                                                                       163225, 98797, 163246,
                                                                       56512, 56542, 100012,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 164737, 0, 3,
                                                                       163246, 98812, 163267,
                                                                       56542, 56572, 100057,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 164800, 0, 3,
                                                                       163288, 98932, 163351,
                                                                       56632, 56692, 100282,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 164926, 0, 3,
                                                                       163351, 98977, 163414,
                                                                       56692, 56752, 100372,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 165052, 0, 3,
                                                                       163414, 99022, 163477,
                                                                       56752, 56812, 100462,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 165178, 0, 3,
                                                                       163477, 99067, 163540,
                                                                       56812, 56872, 100552,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 165304, 0, 3,
                                                                       163540, 99112, 163603,
                                                                       56872, 56932, 100642,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 165430, 0, 3,
                                                                       163603, 99157, 163666,
                                                                       56932, 56992, 100732,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 165556, 0, 3,
                                                                       163666, 99202, 163729,
                                                                       56992, 57052, 100822,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 165682, 0, 3,
                                                                       163729, 99247, 163792,
                                                                       57052, 57112, 100912,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 165808, 0, 3,
                                                                       163792, 99292, 163855,
                                                                       57112, 57172, 101002,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 165934, 0, 3,
                                                                       163855, 99337, 163918,
                                                                       57172, 57232, 101092,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 166060, 0, 3,
                                                                       163918, 99382, 163981,
                                                                       57232, 57292, 101182,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 166186, 0, 3,
                                                                       164044, 99562, 164107,
                                                                       57412, 57472, 101452,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 166312, 0, 3,
                                                                       164107, 99607, 164170,
                                                                       57472, 57532, 101542,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 166438, 0, 3,
                                                                       164170, 99652, 164233,
                                                                       57532, 57592, 101632,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 166564, 0, 3,
                                                                       164233, 99697, 164296,
                                                                       57592, 57652, 101722,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 166690, 0, 3,
                                                                       164296, 99742, 164359,
                                                                       57652, 57712, 101812,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 166816, 0, 3,
                                                                       164359, 99787, 164422,
                                                                       57712, 57772, 101902,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 166942, 0, 3,
                                                                       164422, 99832, 164485,
                                                                       57772, 57832, 101992,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 167068, 0, 3,
                                                                       164485, 99877, 164548,
                                                                       57832, 57892, 102082,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 167194, 0, 3,
                                                                       164548, 99922, 164611,
                                                                       57892, 57952, 102172,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 167320, 0, 3,
                                                                       164611, 99967, 164674,
                                                                       57952, 58012, 102262,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 167446, 0, 3,
                                                                       164674, 100012, 164737,
                                                                       58012, 58072, 102352,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 167572, 0, 3,
                                                                       164800, 100282, 164926,
                                                                       58192, 58292, 102742,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 167782, 0, 3,
                                                                       164926, 100372, 165052,
                                                                       58292, 58392, 102892,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 167992, 0, 3,
                                                                       165052, 100462, 165178,
                                                                       58392, 58492, 103042,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 168202, 0, 3,
                                                                       165178, 100552, 165304,
                                                                       58492, 58592, 103192,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 168412, 0, 3,
                                                                       165304, 100642, 165430,
                                                                       58592, 58692, 103342,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 168622, 0, 3,
                                                                       165430, 100732, 165556,
                                                                       58692, 58792, 103492,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 168832, 0, 3,
                                                                       165556, 100822, 165682,
                                                                       58792, 58892, 103642,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 169042, 0, 3,
                                                                       165682, 100912, 165808,
                                                                       58892, 58992, 103792,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 169252, 0, 3,
                                                                       165808, 101002, 165934,
                                                                       58992, 59092, 103942,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 169462, 0, 3,
                                                                       165934, 101092, 166060,
                                                                       59092, 59192, 104092,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 169672, 0, 3,
                                                                       166186, 101452, 166312,
                                                                       59392, 59492, 104542,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 169882, 0, 3,
                                                                       166312, 101542, 166438,
                                                                       59492, 59592, 104692,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 170092, 0, 3,
                                                                       166438, 101632, 166564,
                                                                       59592, 59692, 104842,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 170302, 0, 3,
                                                                       166564, 101722, 166690,
                                                                       59692, 59792, 104992,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 170512, 0, 3,
                                                                       166690, 101812, 166816,
                                                                       59792, 59892, 105142,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 170722, 0, 3,
                                                                       166816, 101902, 166942,
                                                                       59892, 59992, 105292,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 170932, 0, 3,
                                                                       166942, 101992, 167068,
                                                                       59992, 60092, 105442,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 171142, 0, 3,
                                                                       167068, 102082, 167194,
                                                                       60092, 60192, 105592,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 171352, 0, 3,
                                                                       167194, 102172, 167320,
                                                                       60192, 60292, 105742,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 171562, 0, 3,
                                                                       167320, 102262, 167446,
                                                                       60292, 60392, 105892,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 171772, 0, 3,
                                                                       167572, 102742, 167782,
                                                                       60592, 60742, 106492,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 172087, 0, 3,
                                                                       167782, 102892, 167992,
                                                                       60742, 60892, 106717,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 172402, 0, 3,
                                                                       167992, 103042, 168202,
                                                                       60892, 61042, 106942,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 172717, 0, 3,
                                                                       168202, 103192, 168412,
                                                                       61042, 61192, 107167,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 173032, 0, 3,
                                                                       168412, 103342, 168622,
                                                                       61192, 61342, 107392,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 173347, 0, 3,
                                                                       168622, 103492, 168832,
                                                                       61342, 61492, 107617,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 173662, 0, 3,
                                                                       168832, 103642, 169042,
                                                                       61492, 61642, 107842,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 173977, 0, 3,
                                                                       169042, 103792, 169252,
                                                                       61642, 61792, 108067,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 174292, 0, 3,
                                                                       169252, 103942, 169462,
                                                                       61792, 61942, 108292,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 174607, 0, 3,
                                                                       169672, 104542, 169882,
                                                                       62242, 62392, 108967,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 174922, 0, 3,
                                                                       169882, 104692, 170092,
                                                                       62392, 62542, 109192,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 175237, 0, 3,
                                                                       170092, 104842, 170302,
                                                                       62542, 62692, 109417,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 175552, 0, 3,
                                                                       170302, 104992, 170512,
                                                                       62692, 62842, 109642,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 175867, 0, 3,
                                                                       170512, 105142, 170722,
                                                                       62842, 62992, 109867,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 176182, 0, 3,
                                                                       170722, 105292, 170932,
                                                                       62992, 63142, 110092,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 176497, 0, 3,
                                                                       170932, 105442, 171142,
                                                                       63142, 63292, 110317,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 176812, 0, 3,
                                                                       171142, 105592, 171352,
                                                                       63292, 63442, 110542,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 177127, 0, 3,
                                                                       171352, 105742, 171562,
                                                                       63442, 63592, 110767,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 177442, 0, 3,
                                                                       171772, 106492, 172087,
                                                                       63892, 64102, 111622,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 177883, 0, 3,
                                                                       172087, 106717, 172402,
                                                                       64102, 64312, 111937,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 178324, 0, 3,
                                                                       172402, 106942, 172717,
                                                                       64312, 64522, 112252,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 178765, 0, 3,
                                                                       172717, 107167, 173032,
                                                                       64522, 64732, 112567,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 179206, 0, 3,
                                                                       173032, 107392, 173347,
                                                                       64732, 64942, 112882,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 179647, 0, 3,
                                                                       173347, 107617, 173662,
                                                                       64942, 65152, 113197,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 180088, 0, 3,
                                                                       173662, 107842, 173977,
                                                                       65152, 65362, 113512,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 180529, 0, 3,
                                                                       173977, 108067, 174292,
                                                                       65362, 65572, 113827,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 180970, 0, 3,
                                                                       174607, 108967, 174922,
                                                                       65992, 66202, 114772,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 181411, 0, 3,
                                                                       174922, 109192, 175237,
                                                                       66202, 66412, 115087,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 181852, 0, 3,
                                                                       175237, 109417, 175552,
                                                                       66412, 66622, 115402,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 182293, 0, 3,
                                                                       175552, 109642, 175867,
                                                                       66622, 66832, 115717,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 182734, 0, 3,
                                                                       175867, 109867, 176182,
                                                                       66832, 67042, 116032,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 183175, 0, 3,
                                                                       176182, 110092, 176497,
                                                                       67042, 67252, 116347,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 183616, 0, 3,
                                                                       176497, 110317, 176812,
                                                                       67252, 67462, 116662,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 184057, 0, 3,
                                                                       176812, 110542, 177127,
                                                                       67462, 67672, 116977,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 184498, 0, 3,
                                                                       177442, 111622, 177883,
                                                                       68092, 68372, 118132,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 185086, 0, 3,
                                                                       177883, 111937, 178324,
                                                                       68372, 68652, 118552,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 185674, 0, 3,
                                                                       178324, 112252, 178765,
                                                                       68652, 68932, 118972,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 186262, 0, 3,
                                                                       178765, 112567, 179206,
                                                                       68932, 69212, 119392,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 186850, 0, 3,
                                                                       179206, 112882, 179647,
                                                                       69212, 69492, 119812,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 187438, 0, 3,
                                                                       179647, 113197, 180088,
                                                                       69492, 69772, 120232,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 188026, 0, 3,
                                                                       180088, 113512, 180529,
                                                                       69772, 70052, 120652,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 188614, 0, 3,
                                                                       180970, 114772, 181411,
                                                                       70612, 70892, 121912,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 189202, 0, 3,
                                                                       181411, 115087, 181852,
                                                                       70892, 71172, 122332,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 189790, 0, 3,
                                                                       181852, 115402, 182293,
                                                                       71172, 71452, 122752,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 190378, 0, 3,
                                                                       182293, 115717, 182734,
                                                                       71452, 71732, 123172,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 190966, 0, 3,
                                                                       182734, 116032, 183175,
                                                                       71732, 72012, 123592,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 191554, 0, 3,
                                                                       183175, 116347, 183616,
                                                                       72012, 72292, 124012,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 192142, 0, 3,
                                                                       183616, 116662, 184057,
                                                                       72292, 72572, 124432,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 192730, 0, 3,
                                                                       184498, 118132, 185086,
                                                                       73132, 73492, 125932,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 193486, 0, 3,
                                                                       185086, 118552, 185674,
                                                                       73492, 73852, 126472,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 194242, 0, 3,
                                                                       185674, 118972, 186262,
                                                                       73852, 74212, 127012,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 194998, 0, 3,
                                                                       186262, 119392, 186850,
                                                                       74212, 74572, 127552,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 195754, 0, 3,
                                                                       186850, 119812, 187438,
                                                                       74572, 74932, 128092,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 196510, 0, 3,
                                                                       187438, 120232, 188026,
                                                                       74932, 75292, 128632,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 197266, 0, 3,
                                                                       188614, 121912, 189202,
                                                                       76012, 76372, 130252,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 198022, 0, 3,
                                                                       189202, 122332, 189790,
                                                                       76372, 76732, 130792,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 198778, 0, 3,
                                                                       189790, 122752, 190378,
                                                                       76732, 77092, 131332,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 199534, 0, 3,
                                                                       190378, 123172, 190966,
                                                                       77092, 77452, 131872,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 200290, 0, 3,
                                                                       190966, 123592, 191554,
                                                                       77452, 77812, 132412,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 201046, 0, 3,
                                                                       191554, 124012, 192142,
                                                                       77812, 78172, 132952,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 201802, 0, 3,
                                                                       192730, 125932, 193486,
                                                                       78892, 79342, 134842,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 202747, 0, 3,
                                                                       193486, 126472, 194242,
                                                                       79342, 79792, 135517,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 203692, 0, 3,
                                                                       194242, 127012, 194998,
                                                                       79792, 80242, 136192,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 204637, 0, 3,
                                                                       194998, 127552, 195754,
                                                                       80242, 80692, 136867,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 205582, 0, 3,
                                                                       195754, 128092, 196510,
                                                                       80692, 81142, 137542,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 206527, 0, 3,
                                                                       197266, 130252, 198022,
                                                                       82042, 82492, 139567,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 207472, 0, 3,
                                                                       198022, 130792, 198778,
                                                                       82492, 82942, 140242,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 208417, 0, 3,
                                                                       198778, 131332, 199534,
                                                                       82942, 83392, 140917,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 209362, 0, 3,
                                                                       199534, 131872, 200290,
                                                                       83392, 83842, 141592,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 210307, 0, 3,
                                                                       200290, 132412, 201046,
                                                                       83842, 84292, 142267,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 211252, 0, 3,
                                                                       201802, 134842, 202747,
                                                                       85192, 85742, 144592,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 212407, 0, 3,
                                                                       202747, 135517, 203692,
                                                                       85742, 86292, 145417,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 213562, 0, 3,
                                                                       203692, 136192, 204637,
                                                                       86292, 86842, 146242,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 214717, 0, 3,
                                                                       204637, 136867, 205582,
                                                                       86842, 87392, 147067,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 215872, 0, 3,
                                                                       206527, 139567, 207472,
                                                                       88492, 89042, 149542,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 217027, 0, 3,
                                                                       207472, 140242, 208417,
                                                                       89042, 89592, 150367,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 218182, 0, 3,
                                                                       208417, 140917, 209362,
                                                                       89592, 90142, 151192,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 219337, 0, 3,
                                                                       209362, 141592, 210307,
                                                                       90142, 90692, 152017,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 220492, 0, 3,
                                                                       211252, 144592, 212407,
                                                                       91792, 92452, 154822,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 221878, 0, 3,
                                                                       212407, 145417, 213562,
                                                                       92452, 93112, 155812,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 223264, 0, 3,
                                                                       213562, 146242, 214717,
                                                                       93112, 93772, 156802,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 224650, 0, 3,
                                                                       215872, 149542, 217027,
                                                                       95092, 95752, 159772,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 226036, 0, 3,
                                                                       217027, 150367, 218182,
                                                                       95752, 96412, 160762,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 227422, 0, 3,
                                                                       218182, 151192, 219337,
                                                                       96412, 97072, 161752,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 228808, 3, 98392,
                                                                       98407, 162742, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 228836, 3, 98407,
                                                                       98422, 162763, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 228864, 3, 98422,
                                                                       98437, 162784, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 228892, 3, 98437,
                                                                       98452, 162805, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 228920, 3, 98452,
                                                                       98467, 162826, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 228948, 3, 98467,
                                                                       98482, 162847, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 228976, 3, 98482,
                                                                       98497, 162868, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 229004, 3, 98497,
                                                                       98512, 162889, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 229032, 3, 98512,
                                                                       98527, 162910, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 229060, 3, 98527,
                                                                       98542, 162931, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 229088, 3, 98542,
                                                                       98557, 162952, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 229116, 3, 98557,
                                                                       98572, 162973, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 229144, 3, 98572,
                                                                       98587, 162994, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 229172, 3, 98617,
                                                                       98632, 163015, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 229200, 3, 98632,
                                                                       98647, 163036, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 229228, 3, 98647,
                                                                       98662, 163057, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 229256, 3, 98662,
                                                                       98677, 163078, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 229284, 3, 98677,
                                                                       98692, 163099, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 229312, 3, 98692,
                                                                       98707, 163120, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 229340, 3, 98707,
                                                                       98722, 163141, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 229368, 3, 98722,
                                                                       98737, 163162, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 229396, 3, 98737,
                                                                       98752, 163183, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 229424, 3, 98752,
                                                                       98767, 163204, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 229452, 3, 98767,
                                                                       98782, 163225, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 229480, 3, 98782,
                                                                       98797, 163246, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 229508, 3, 98797,
                                                                       98812, 163267, ncols,
                                                                       gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 229536, 0, 3,
                                                                       228808, 162742, 228836,
                                                                       98842, 98887, 163288,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 229620, 0, 3,
                                                                       228836, 162763, 228864,
                                                                       98887, 98932, 163351,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 229704, 0, 3,
                                                                       228864, 162784, 228892,
                                                                       98932, 98977, 163414,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 229788, 0, 3,
                                                                       228892, 162805, 228920,
                                                                       98977, 99022, 163477,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 229872, 0, 3,
                                                                       228920, 162826, 228948,
                                                                       99022, 99067, 163540,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 229956, 0, 3,
                                                                       228948, 162847, 228976,
                                                                       99067, 99112, 163603,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 230040, 0, 3,
                                                                       228976, 162868, 229004,
                                                                       99112, 99157, 163666,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 230124, 0, 3,
                                                                       229004, 162889, 229032,
                                                                       99157, 99202, 163729,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 230208, 0, 3,
                                                                       229032, 162910, 229060,
                                                                       99202, 99247, 163792,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 230292, 0, 3,
                                                                       229060, 162931, 229088,
                                                                       99247, 99292, 163855,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 230376, 0, 3,
                                                                       229088, 162952, 229116,
                                                                       99292, 99337, 163918,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 230460, 0, 3,
                                                                       229116, 162973, 229144,
                                                                       99337, 99382, 163981,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 230544, 0, 3,
                                                                       229172, 163015, 229200,
                                                                       99472, 99517, 164044,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 230628, 0, 3,
                                                                       229200, 163036, 229228,
                                                                       99517, 99562, 164107,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 230712, 0, 3,
                                                                       229228, 163057, 229256,
                                                                       99562, 99607, 164170,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 230796, 0, 3,
                                                                       229256, 163078, 229284,
                                                                       99607, 99652, 164233,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 230880, 0, 3,
                                                                       229284, 163099, 229312,
                                                                       99652, 99697, 164296,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 230964, 0, 3,
                                                                       229312, 163120, 229340,
                                                                       99697, 99742, 164359,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 231048, 0, 3,
                                                                       229340, 163141, 229368,
                                                                       99742, 99787, 164422,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 231132, 0, 3,
                                                                       229368, 163162, 229396,
                                                                       99787, 99832, 164485,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 231216, 0, 3,
                                                                       229396, 163183, 229424,
                                                                       99832, 99877, 164548,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 231300, 0, 3,
                                                                       229424, 163204, 229452,
                                                                       99877, 99922, 164611,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 231384, 0, 3,
                                                                       229452, 163225, 229480,
                                                                       99922, 99967, 164674,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 231468, 0, 3,
                                                                       229480, 163246, 229508,
                                                                       99967, 100012, 164737,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 231552, 0, 3,
                                                                       229536, 163288, 229620,
                                                                       100102, 100192, 164800,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 231720, 0, 3,
                                                                       229620, 163351, 229704,
                                                                       100192, 100282, 164926,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 231888, 0, 3,
                                                                       229704, 163414, 229788,
                                                                       100282, 100372, 165052,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 232056, 0, 3,
                                                                       229788, 163477, 229872,
                                                                       100372, 100462, 165178,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 232224, 0, 3,
                                                                       229872, 163540, 229956,
                                                                       100462, 100552, 165304,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 232392, 0, 3,
                                                                       229956, 163603, 230040,
                                                                       100552, 100642, 165430,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 232560, 0, 3,
                                                                       230040, 163666, 230124,
                                                                       100642, 100732, 165556,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 232728, 0, 3,
                                                                       230124, 163729, 230208,
                                                                       100732, 100822, 165682,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 232896, 0, 3,
                                                                       230208, 163792, 230292,
                                                                       100822, 100912, 165808,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 233064, 0, 3,
                                                                       230292, 163855, 230376,
                                                                       100912, 101002, 165934,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 233232, 0, 3,
                                                                       230376, 163918, 230460,
                                                                       101002, 101092, 166060,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 233400, 0, 3,
                                                                       230544, 164044, 230628,
                                                                       101272, 101362, 166186,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 233568, 0, 3,
                                                                       230628, 164107, 230712,
                                                                       101362, 101452, 166312,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 233736, 0, 3,
                                                                       230712, 164170, 230796,
                                                                       101452, 101542, 166438,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 233904, 0, 3,
                                                                       230796, 164233, 230880,
                                                                       101542, 101632, 166564,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 234072, 0, 3,
                                                                       230880, 164296, 230964,
                                                                       101632, 101722, 166690,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 234240, 0, 3,
                                                                       230964, 164359, 231048,
                                                                       101722, 101812, 166816,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 234408, 0, 3,
                                                                       231048, 164422, 231132,
                                                                       101812, 101902, 166942,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 234576, 0, 3,
                                                                       231132, 164485, 231216,
                                                                       101902, 101992, 167068,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 234744, 0, 3,
                                                                       231216, 164548, 231300,
                                                                       101992, 102082, 167194,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 234912, 0, 3,
                                                                       231300, 164611, 231384,
                                                                       102082, 102172, 167320,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 235080, 0, 3,
                                                                       231384, 164674, 231468,
                                                                       102172, 102262, 167446,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 235248, 0, 3,
                                                                       231552, 164800, 231720,
                                                                       102442, 102592, 167572,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 235528, 0, 3,
                                                                       231720, 164926, 231888,
                                                                       102592, 102742, 167782,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 235808, 0, 3,
                                                                       231888, 165052, 232056,
                                                                       102742, 102892, 167992,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 236088, 0, 3,
                                                                       232056, 165178, 232224,
                                                                       102892, 103042, 168202,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 236368, 0, 3,
                                                                       232224, 165304, 232392,
                                                                       103042, 103192, 168412,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 236648, 0, 3,
                                                                       232392, 165430, 232560,
                                                                       103192, 103342, 168622,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 236928, 0, 3,
                                                                       232560, 165556, 232728,
                                                                       103342, 103492, 168832,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 237208, 0, 3,
                                                                       232728, 165682, 232896,
                                                                       103492, 103642, 169042,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 237488, 0, 3,
                                                                       232896, 165808, 233064,
                                                                       103642, 103792, 169252,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 237768, 0, 3,
                                                                       233064, 165934, 233232,
                                                                       103792, 103942, 169462,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 238048, 0, 3,
                                                                       233400, 166186, 233568,
                                                                       104242, 104392, 169672,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 238328, 0, 3,
                                                                       233568, 166312, 233736,
                                                                       104392, 104542, 169882,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 238608, 0, 3,
                                                                       233736, 166438, 233904,
                                                                       104542, 104692, 170092,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 238888, 0, 3,
                                                                       233904, 166564, 234072,
                                                                       104692, 104842, 170302,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 239168, 0, 3,
                                                                       234072, 166690, 234240,
                                                                       104842, 104992, 170512,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 239448, 0, 3,
                                                                       234240, 166816, 234408,
                                                                       104992, 105142, 170722,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 239728, 0, 3,
                                                                       234408, 166942, 234576,
                                                                       105142, 105292, 170932,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 240008, 0, 3,
                                                                       234576, 167068, 234744,
                                                                       105292, 105442, 171142,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 240288, 0, 3,
                                                                       234744, 167194, 234912,
                                                                       105442, 105592, 171352,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 240568, 0, 3,
                                                                       234912, 167320, 235080,
                                                                       105592, 105742, 171562,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 240848, 0, 3,
                                                                       235248, 167572, 235528,
                                                                       106042, 106267, 171772,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 241268, 0, 3,
                                                                       235528, 167782, 235808,
                                                                       106267, 106492, 172087,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 241688, 0, 3,
                                                                       235808, 167992, 236088,
                                                                       106492, 106717, 172402,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 242108, 0, 3,
                                                                       236088, 168202, 236368,
                                                                       106717, 106942, 172717,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 242528, 0, 3,
                                                                       236368, 168412, 236648,
                                                                       106942, 107167, 173032,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 242948, 0, 3,
                                                                       236648, 168622, 236928,
                                                                       107167, 107392, 173347,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 243368, 0, 3,
                                                                       236928, 168832, 237208,
                                                                       107392, 107617, 173662,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 243788, 0, 3,
                                                                       237208, 169042, 237488,
                                                                       107617, 107842, 173977,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 244208, 0, 3,
                                                                       237488, 169252, 237768,
                                                                       107842, 108067, 174292,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 244628, 0, 3,
                                                                       238048, 169672, 238328,
                                                                       108517, 108742, 174607,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 245048, 0, 3,
                                                                       238328, 169882, 238608,
                                                                       108742, 108967, 174922,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 245468, 0, 3,
                                                                       238608, 170092, 238888,
                                                                       108967, 109192, 175237,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 245888, 0, 3,
                                                                       238888, 170302, 239168,
                                                                       109192, 109417, 175552,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 246308, 0, 3,
                                                                       239168, 170512, 239448,
                                                                       109417, 109642, 175867,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 246728, 0, 3,
                                                                       239448, 170722, 239728,
                                                                       109642, 109867, 176182,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 247148, 0, 3,
                                                                       239728, 170932, 240008,
                                                                       109867, 110092, 176497,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 247568, 0, 3,
                                                                       240008, 171142, 240288,
                                                                       110092, 110317, 176812,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 247988, 0, 3,
                                                                       240288, 171352, 240568,
                                                                       110317, 110542, 177127,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 248408, 0, 3,
                                                                       240848, 171772, 241268,
                                                                       110992, 111307, 177442,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 248996, 0, 3,
                                                                       241268, 172087, 241688,
                                                                       111307, 111622, 177883,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 249584, 0, 3,
                                                                       241688, 172402, 242108,
                                                                       111622, 111937, 178324,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 250172, 0, 3,
                                                                       242108, 172717, 242528,
                                                                       111937, 112252, 178765,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 250760, 0, 3,
                                                                       242528, 173032, 242948,
                                                                       112252, 112567, 179206,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 251348, 0, 3,
                                                                       242948, 173347, 243368,
                                                                       112567, 112882, 179647,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 251936, 0, 3,
                                                                       243368, 173662, 243788,
                                                                       112882, 113197, 180088,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 252524, 0, 3,
                                                                       243788, 173977, 244208,
                                                                       113197, 113512, 180529,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 253112, 0, 3,
                                                                       244628, 174607, 245048,
                                                                       114142, 114457, 180970,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 253700, 0, 3,
                                                                       245048, 174922, 245468,
                                                                       114457, 114772, 181411,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 254288, 0, 3,
                                                                       245468, 175237, 245888,
                                                                       114772, 115087, 181852,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 254876, 0, 3,
                                                                       245888, 175552, 246308,
                                                                       115087, 115402, 182293,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 255464, 0, 3,
                                                                       246308, 175867, 246728,
                                                                       115402, 115717, 182734,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 256052, 0, 3,
                                                                       246728, 176182, 247148,
                                                                       115717, 116032, 183175,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 256640, 0, 3,
                                                                       247148, 176497, 247568,
                                                                       116032, 116347, 183616,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 257228, 0, 3,
                                                                       247568, 176812, 247988,
                                                                       116347, 116662, 184057,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 257816, 0, 3,
                                                                       248408, 177442, 248996,
                                                                       117292, 117712, 184498,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 258600, 0, 3,
                                                                       248996, 177883, 249584,
                                                                       117712, 118132, 185086,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 259384, 0, 3,
                                                                       249584, 178324, 250172,
                                                                       118132, 118552, 185674,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 260168, 0, 3,
                                                                       250172, 178765, 250760,
                                                                       118552, 118972, 186262,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 260952, 0, 3,
                                                                       250760, 179206, 251348,
                                                                       118972, 119392, 186850,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 261736, 0, 3,
                                                                       251348, 179647, 251936,
                                                                       119392, 119812, 187438,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 262520, 0, 3,
                                                                       251936, 180088, 252524,
                                                                       119812, 120232, 188026,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 263304, 0, 3,
                                                                       253112, 180970, 253700,
                                                                       121072, 121492, 188614,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 264088, 0, 3,
                                                                       253700, 181411, 254288,
                                                                       121492, 121912, 189202,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 264872, 0, 3,
                                                                       254288, 181852, 254876,
                                                                       121912, 122332, 189790,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 265656, 0, 3,
                                                                       254876, 182293, 255464,
                                                                       122332, 122752, 190378,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 266440, 0, 3,
                                                                       255464, 182734, 256052,
                                                                       122752, 123172, 190966,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 267224, 0, 3,
                                                                       256052, 183175, 256640,
                                                                       123172, 123592, 191554,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 268008, 0, 3,
                                                                       256640, 183616, 257228,
                                                                       123592, 124012, 192142,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 268792, 0, 3,
                                                                       257816, 184498, 258600,
                                                                       124852, 125392, 192730,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 269800, 0, 3,
                                                                       258600, 185086, 259384,
                                                                       125392, 125932, 193486,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 270808, 0, 3,
                                                                       259384, 185674, 260168,
                                                                       125932, 126472, 194242,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 271816, 0, 3,
                                                                       260168, 186262, 260952,
                                                                       126472, 127012, 194998,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 272824, 0, 3,
                                                                       260952, 186850, 261736,
                                                                       127012, 127552, 195754,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 273832, 0, 3,
                                                                       261736, 187438, 262520,
                                                                       127552, 128092, 196510,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 274840, 0, 3,
                                                                       263304, 188614, 264088,
                                                                       129172, 129712, 197266,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 275848, 0, 3,
                                                                       264088, 189202, 264872,
                                                                       129712, 130252, 198022,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 276856, 0, 3,
                                                                       264872, 189790, 265656,
                                                                       130252, 130792, 198778,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 277864, 0, 3,
                                                                       265656, 190378, 266440,
                                                                       130792, 131332, 199534,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 278872, 0, 3,
                                                                       266440, 190966, 267224,
                                                                       131332, 131872, 200290,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 279880, 0, 3,
                                                                       267224, 191554, 268008,
                                                                       131872, 132412, 201046,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 280888, 0, 3,
                                                                       268792, 192730, 269800,
                                                                       133492, 134167, 201802,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 282148, 0, 3,
                                                                       269800, 193486, 270808,
                                                                       134167, 134842, 202747,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 283408, 0, 3,
                                                                       270808, 194242, 271816,
                                                                       134842, 135517, 203692,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 284668, 0, 3,
                                                                       271816, 194998, 272824,
                                                                       135517, 136192, 204637,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 285928, 0, 3,
                                                                       272824, 195754, 273832,
                                                                       136192, 136867, 205582,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 287188, 0, 3,
                                                                       274840, 197266, 275848,
                                                                       138217, 138892, 206527,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 288448, 0, 3,
                                                                       275848, 198022, 276856,
                                                                       138892, 139567, 207472,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 289708, 0, 3,
                                                                       276856, 198778, 277864,
                                                                       139567, 140242, 208417,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 290968, 0, 3,
                                                                       277864, 199534, 278872,
                                                                       140242, 140917, 209362,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 292228, 0, 3,
                                                                       278872, 200290, 279880,
                                                                       140917, 141592, 210307,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 293488, 0, 3,
                                                                       280888, 201802, 282148,
                                                                       142942, 143767, 211252,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 295028, 0, 3,
                                                                       282148, 202747, 283408,
                                                                       143767, 144592, 212407,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 296568, 0, 3,
                                                                       283408, 203692, 284668,
                                                                       144592, 145417, 213562,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 298108, 0, 3,
                                                                       284668, 204637, 285928,
                                                                       145417, 146242, 214717,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 299648, 0, 3,
                                                                       287188, 206527, 288448,
                                                                       147892, 148717, 215872,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 301188, 0, 3,
                                                                       288448, 207472, 289708,
                                                                       148717, 149542, 217027,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 302728, 0, 3,
                                                                       289708, 208417, 290968,
                                                                       149542, 150367, 218182,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 304268, 0, 3,
                                                                       290968, 209362, 292228,
                                                                       150367, 151192, 219337,
                                                                       ncols, gamma, p, q);

                    compute_prim_sni_three_center_electron_repulsion_0(buffer, 305808, 0, 3,
                                                                       293488, 211252, 295028,
                                                                       152842, 153832, 220492,
                                                                       ncols, gamma, p, q);

                    compute_prim_sni_three_center_electron_repulsion_0(buffer, 307656, 0, 3,
                                                                       295028, 212407, 296568,
                                                                       153832, 154822, 221878,
                                                                       ncols, gamma, p, q);

                    compute_prim_sni_three_center_electron_repulsion_0(buffer, 309504, 0, 3,
                                                                       296568, 213562, 298108,
                                                                       154822, 155812, 223264,
                                                                       ncols, gamma, p, q);

                    compute_prim_sni_three_center_electron_repulsion_0(buffer, 311352, 0, 3,
                                                                       299648, 215872, 301188,
                                                                       157792, 158782, 224650,
                                                                       ncols, gamma, p, q);

                    compute_prim_sni_three_center_electron_repulsion_0(buffer, 313200, 0, 3,
                                                                       301188, 217027, 302728,
                                                                       158782, 159772, 226036,
                                                                       ncols, gamma, p, q);

                    compute_prim_sni_three_center_electron_repulsion_0(buffer, 315048, 0, 3,
                                                                       302728, 218182, 304268,
                                                                       159772, 160762, 227422,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 316896, 3, 162742,
                                                                       162763, 228864, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 316932, 3, 162763,
                                                                       162784, 228892, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 316968, 3, 162784,
                                                                       162805, 228920, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 317004, 3, 162805,
                                                                       162826, 228948, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 317040, 3, 162826,
                                                                       162847, 228976, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 317076, 3, 162847,
                                                                       162868, 229004, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 317112, 3, 162868,
                                                                       162889, 229032, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 317148, 3, 162889,
                                                                       162910, 229060, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 317184, 3, 162910,
                                                                       162931, 229088, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 317220, 3, 162931,
                                                                       162952, 229116, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 317256, 3, 162952,
                                                                       162973, 229144, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 317292, 3, 163015,
                                                                       163036, 229228, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 317328, 3, 163036,
                                                                       163057, 229256, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 317364, 3, 163057,
                                                                       163078, 229284, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 317400, 3, 163078,
                                                                       163099, 229312, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 317436, 3, 163099,
                                                                       163120, 229340, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 317472, 3, 163120,
                                                                       163141, 229368, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 317508, 3, 163141,
                                                                       163162, 229396, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 317544, 3, 163162,
                                                                       163183, 229424, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 317580, 3, 163183,
                                                                       163204, 229452, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 317616, 3, 163204,
                                                                       163225, 229480, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 317652, 3, 163225,
                                                                       163246, 229508, ncols,
                                                                       gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 317688, 0, 3,
                                                                       316896, 228864, 316932,
                                                                       163288, 163351, 229704,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 317796, 0, 3,
                                                                       316932, 228892, 316968,
                                                                       163351, 163414, 229788,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 317904, 0, 3,
                                                                       316968, 228920, 317004,
                                                                       163414, 163477, 229872,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 318012, 0, 3,
                                                                       317004, 228948, 317040,
                                                                       163477, 163540, 229956,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 318120, 0, 3,
                                                                       317040, 228976, 317076,
                                                                       163540, 163603, 230040,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 318228, 0, 3,
                                                                       317076, 229004, 317112,
                                                                       163603, 163666, 230124,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 318336, 0, 3,
                                                                       317112, 229032, 317148,
                                                                       163666, 163729, 230208,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 318444, 0, 3,
                                                                       317148, 229060, 317184,
                                                                       163729, 163792, 230292,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 318552, 0, 3,
                                                                       317184, 229088, 317220,
                                                                       163792, 163855, 230376,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 318660, 0, 3,
                                                                       317220, 229116, 317256,
                                                                       163855, 163918, 230460,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 318768, 0, 3,
                                                                       317292, 229228, 317328,
                                                                       164044, 164107, 230712,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 318876, 0, 3,
                                                                       317328, 229256, 317364,
                                                                       164107, 164170, 230796,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 318984, 0, 3,
                                                                       317364, 229284, 317400,
                                                                       164170, 164233, 230880,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 319092, 0, 3,
                                                                       317400, 229312, 317436,
                                                                       164233, 164296, 230964,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 319200, 0, 3,
                                                                       317436, 229340, 317472,
                                                                       164296, 164359, 231048,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 319308, 0, 3,
                                                                       317472, 229368, 317508,
                                                                       164359, 164422, 231132,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 319416, 0, 3,
                                                                       317508, 229396, 317544,
                                                                       164422, 164485, 231216,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 319524, 0, 3,
                                                                       317544, 229424, 317580,
                                                                       164485, 164548, 231300,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 319632, 0, 3,
                                                                       317580, 229452, 317616,
                                                                       164548, 164611, 231384,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 319740, 0, 3,
                                                                       317616, 229480, 317652,
                                                                       164611, 164674, 231468,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 319848, 0, 3,
                                                                       317688, 229704, 317796,
                                                                       164800, 164926, 231888,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 320064, 0, 3,
                                                                       317796, 229788, 317904,
                                                                       164926, 165052, 232056,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 320280, 0, 3,
                                                                       317904, 229872, 318012,
                                                                       165052, 165178, 232224,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 320496, 0, 3,
                                                                       318012, 229956, 318120,
                                                                       165178, 165304, 232392,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 320712, 0, 3,
                                                                       318120, 230040, 318228,
                                                                       165304, 165430, 232560,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 320928, 0, 3,
                                                                       318228, 230124, 318336,
                                                                       165430, 165556, 232728,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 321144, 0, 3,
                                                                       318336, 230208, 318444,
                                                                       165556, 165682, 232896,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 321360, 0, 3,
                                                                       318444, 230292, 318552,
                                                                       165682, 165808, 233064,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 321576, 0, 3,
                                                                       318552, 230376, 318660,
                                                                       165808, 165934, 233232,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 321792, 0, 3,
                                                                       318768, 230712, 318876,
                                                                       166186, 166312, 233736,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 322008, 0, 3,
                                                                       318876, 230796, 318984,
                                                                       166312, 166438, 233904,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 322224, 0, 3,
                                                                       318984, 230880, 319092,
                                                                       166438, 166564, 234072,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 322440, 0, 3,
                                                                       319092, 230964, 319200,
                                                                       166564, 166690, 234240,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 322656, 0, 3,
                                                                       319200, 231048, 319308,
                                                                       166690, 166816, 234408,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 322872, 0, 3,
                                                                       319308, 231132, 319416,
                                                                       166816, 166942, 234576,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 323088, 0, 3,
                                                                       319416, 231216, 319524,
                                                                       166942, 167068, 234744,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 323304, 0, 3,
                                                                       319524, 231300, 319632,
                                                                       167068, 167194, 234912,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 323520, 0, 3,
                                                                       319632, 231384, 319740,
                                                                       167194, 167320, 235080,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 323736, 0, 3,
                                                                       319848, 231888, 320064,
                                                                       167572, 167782, 235808,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 324096, 0, 3,
                                                                       320064, 232056, 320280,
                                                                       167782, 167992, 236088,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 324456, 0, 3,
                                                                       320280, 232224, 320496,
                                                                       167992, 168202, 236368,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 324816, 0, 3,
                                                                       320496, 232392, 320712,
                                                                       168202, 168412, 236648,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 325176, 0, 3,
                                                                       320712, 232560, 320928,
                                                                       168412, 168622, 236928,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 325536, 0, 3,
                                                                       320928, 232728, 321144,
                                                                       168622, 168832, 237208,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 325896, 0, 3,
                                                                       321144, 232896, 321360,
                                                                       168832, 169042, 237488,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 326256, 0, 3,
                                                                       321360, 233064, 321576,
                                                                       169042, 169252, 237768,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 326616, 0, 3,
                                                                       321792, 233736, 322008,
                                                                       169672, 169882, 238608,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 326976, 0, 3,
                                                                       322008, 233904, 322224,
                                                                       169882, 170092, 238888,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 327336, 0, 3,
                                                                       322224, 234072, 322440,
                                                                       170092, 170302, 239168,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 327696, 0, 3,
                                                                       322440, 234240, 322656,
                                                                       170302, 170512, 239448,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 328056, 0, 3,
                                                                       322656, 234408, 322872,
                                                                       170512, 170722, 239728,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 328416, 0, 3,
                                                                       322872, 234576, 323088,
                                                                       170722, 170932, 240008,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 328776, 0, 3,
                                                                       323088, 234744, 323304,
                                                                       170932, 171142, 240288,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 329136, 0, 3,
                                                                       323304, 234912, 323520,
                                                                       171142, 171352, 240568,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 329496, 0, 3,
                                                                       323736, 235808, 324096,
                                                                       171772, 172087, 241688,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 330036, 0, 3,
                                                                       324096, 236088, 324456,
                                                                       172087, 172402, 242108,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 330576, 0, 3,
                                                                       324456, 236368, 324816,
                                                                       172402, 172717, 242528,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 331116, 0, 3,
                                                                       324816, 236648, 325176,
                                                                       172717, 173032, 242948,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 331656, 0, 3,
                                                                       325176, 236928, 325536,
                                                                       173032, 173347, 243368,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 332196, 0, 3,
                                                                       325536, 237208, 325896,
                                                                       173347, 173662, 243788,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 332736, 0, 3,
                                                                       325896, 237488, 326256,
                                                                       173662, 173977, 244208,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 333276, 0, 3,
                                                                       326616, 238608, 326976,
                                                                       174607, 174922, 245468,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 333816, 0, 3,
                                                                       326976, 238888, 327336,
                                                                       174922, 175237, 245888,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 334356, 0, 3,
                                                                       327336, 239168, 327696,
                                                                       175237, 175552, 246308,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 334896, 0, 3,
                                                                       327696, 239448, 328056,
                                                                       175552, 175867, 246728,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 335436, 0, 3,
                                                                       328056, 239728, 328416,
                                                                       175867, 176182, 247148,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 335976, 0, 3,
                                                                       328416, 240008, 328776,
                                                                       176182, 176497, 247568,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 336516, 0, 3,
                                                                       328776, 240288, 329136,
                                                                       176497, 176812, 247988,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 337056, 0, 3,
                                                                       329496, 241688, 330036,
                                                                       177442, 177883, 249584,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 337812, 0, 3,
                                                                       330036, 242108, 330576,
                                                                       177883, 178324, 250172,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 338568, 0, 3,
                                                                       330576, 242528, 331116,
                                                                       178324, 178765, 250760,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 339324, 0, 3,
                                                                       331116, 242948, 331656,
                                                                       178765, 179206, 251348,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 340080, 0, 3,
                                                                       331656, 243368, 332196,
                                                                       179206, 179647, 251936,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 340836, 0, 3,
                                                                       332196, 243788, 332736,
                                                                       179647, 180088, 252524,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 341592, 0, 3,
                                                                       333276, 245468, 333816,
                                                                       180970, 181411, 254288,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 342348, 0, 3,
                                                                       333816, 245888, 334356,
                                                                       181411, 181852, 254876,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 343104, 0, 3,
                                                                       334356, 246308, 334896,
                                                                       181852, 182293, 255464,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 343860, 0, 3,
                                                                       334896, 246728, 335436,
                                                                       182293, 182734, 256052,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 344616, 0, 3,
                                                                       335436, 247148, 335976,
                                                                       182734, 183175, 256640,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 345372, 0, 3,
                                                                       335976, 247568, 336516,
                                                                       183175, 183616, 257228,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 346128, 0, 3,
                                                                       337056, 249584, 337812,
                                                                       184498, 185086, 259384,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 347136, 0, 3,
                                                                       337812, 250172, 338568,
                                                                       185086, 185674, 260168,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 348144, 0, 3,
                                                                       338568, 250760, 339324,
                                                                       185674, 186262, 260952,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 349152, 0, 3,
                                                                       339324, 251348, 340080,
                                                                       186262, 186850, 261736,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 350160, 0, 3,
                                                                       340080, 251936, 340836,
                                                                       186850, 187438, 262520,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 351168, 0, 3,
                                                                       341592, 254288, 342348,
                                                                       188614, 189202, 264872,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 352176, 0, 3,
                                                                       342348, 254876, 343104,
                                                                       189202, 189790, 265656,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 353184, 0, 3,
                                                                       343104, 255464, 343860,
                                                                       189790, 190378, 266440,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 354192, 0, 3,
                                                                       343860, 256052, 344616,
                                                                       190378, 190966, 267224,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 355200, 0, 3,
                                                                       344616, 256640, 345372,
                                                                       190966, 191554, 268008,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 356208, 0, 3,
                                                                       346128, 259384, 347136,
                                                                       192730, 193486, 270808,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 357504, 0, 3,
                                                                       347136, 260168, 348144,
                                                                       193486, 194242, 271816,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 358800, 0, 3,
                                                                       348144, 260952, 349152,
                                                                       194242, 194998, 272824,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 360096, 0, 3,
                                                                       349152, 261736, 350160,
                                                                       194998, 195754, 273832,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 361392, 0, 3,
                                                                       351168, 264872, 352176,
                                                                       197266, 198022, 276856,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 362688, 0, 3,
                                                                       352176, 265656, 353184,
                                                                       198022, 198778, 277864,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 363984, 0, 3,
                                                                       353184, 266440, 354192,
                                                                       198778, 199534, 278872,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 365280, 0, 3,
                                                                       354192, 267224, 355200,
                                                                       199534, 200290, 279880,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 366576, 0, 3,
                                                                       356208, 270808, 357504,
                                                                       201802, 202747, 283408,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 368196, 0, 3,
                                                                       357504, 271816, 358800,
                                                                       202747, 203692, 284668,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 369816, 0, 3,
                                                                       358800, 272824, 360096,
                                                                       203692, 204637, 285928,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 371436, 0, 3,
                                                                       361392, 276856, 362688,
                                                                       206527, 207472, 289708,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 373056, 0, 3,
                                                                       362688, 277864, 363984,
                                                                       207472, 208417, 290968,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 374676, 0, 3,
                                                                       363984, 278872, 365280,
                                                                       208417, 209362, 292228,
                                                                       ncols, gamma, p, q);

                    compute_prim_smk_three_center_electron_repulsion_0(buffer, 376296, 0, 3,
                                                                       366576, 283408, 368196,
                                                                       211252, 212407, 296568,
                                                                       ncols, gamma, p, q);

                    compute_prim_smk_three_center_electron_repulsion_0(buffer, 378276, 0, 3,
                                                                       368196, 284668, 369816,
                                                                       212407, 213562, 298108,
                                                                       ncols, gamma, p, q);

                    compute_prim_smk_three_center_electron_repulsion_0(buffer, 380256, 0, 3,
                                                                       371436, 289708, 373056,
                                                                       215872, 217027, 302728,
                                                                       ncols, gamma, p, q);

                    compute_prim_smk_three_center_electron_repulsion_0(buffer, 382236, 0, 3,
                                                                       373056, 290968, 374676,
                                                                       217027, 218182, 304268,
                                                                       ncols, gamma, p, q);

                    compute_prim_snk_three_center_electron_repulsion_0(buffer, 384216, 0, 3,
                                                                       376296, 296568, 378276,
                                                                       220492, 221878, 309504,
                                                                       ncols, gamma, p, q);

                    compute_prim_snk_three_center_electron_repulsion_0(buffer, 386592, 0, 3,
                                                                       380256, 302728, 382236,
                                                                       224650, 226036, 315048,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 388968, 3, 228808,
                                                                       228836, 316896, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 389013, 3, 228836,
                                                                       228864, 316932, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 389058, 3, 228864,
                                                                       228892, 316968, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 389103, 3, 228892,
                                                                       228920, 317004, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 389148, 3, 228920,
                                                                       228948, 317040, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 389193, 3, 228948,
                                                                       228976, 317076, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 389238, 3, 228976,
                                                                       229004, 317112, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 389283, 3, 229004,
                                                                       229032, 317148, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 389328, 3, 229032,
                                                                       229060, 317184, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 389373, 3, 229060,
                                                                       229088, 317220, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 389418, 3, 229088,
                                                                       229116, 317256, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 389463, 3, 229172,
                                                                       229200, 317292, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 389508, 3, 229200,
                                                                       229228, 317328, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 389553, 3, 229228,
                                                                       229256, 317364, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 389598, 3, 229256,
                                                                       229284, 317400, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 389643, 3, 229284,
                                                                       229312, 317436, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 389688, 3, 229312,
                                                                       229340, 317472, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 389733, 3, 229340,
                                                                       229368, 317508, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 389778, 3, 229368,
                                                                       229396, 317544, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 389823, 3, 229396,
                                                                       229424, 317580, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 389868, 3, 229424,
                                                                       229452, 317616, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 389913, 3, 229452,
                                                                       229480, 317652, ncols,
                                                                       gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 389958, 0, 3,
                                                                       388968, 316896, 389013,
                                                                       229536, 229620, 317688,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 390093, 0, 3,
                                                                       389013, 316932, 389058,
                                                                       229620, 229704, 317796,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 390228, 0, 3,
                                                                       389058, 316968, 389103,
                                                                       229704, 229788, 317904,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 390363, 0, 3,
                                                                       389103, 317004, 389148,
                                                                       229788, 229872, 318012,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 390498, 0, 3,
                                                                       389148, 317040, 389193,
                                                                       229872, 229956, 318120,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 390633, 0, 3,
                                                                       389193, 317076, 389238,
                                                                       229956, 230040, 318228,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 390768, 0, 3,
                                                                       389238, 317112, 389283,
                                                                       230040, 230124, 318336,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 390903, 0, 3,
                                                                       389283, 317148, 389328,
                                                                       230124, 230208, 318444,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 391038, 0, 3,
                                                                       389328, 317184, 389373,
                                                                       230208, 230292, 318552,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 391173, 0, 3,
                                                                       389373, 317220, 389418,
                                                                       230292, 230376, 318660,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 391308, 0, 3,
                                                                       389463, 317292, 389508,
                                                                       230544, 230628, 318768,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 391443, 0, 3,
                                                                       389508, 317328, 389553,
                                                                       230628, 230712, 318876,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 391578, 0, 3,
                                                                       389553, 317364, 389598,
                                                                       230712, 230796, 318984,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 391713, 0, 3,
                                                                       389598, 317400, 389643,
                                                                       230796, 230880, 319092,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 391848, 0, 3,
                                                                       389643, 317436, 389688,
                                                                       230880, 230964, 319200,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 391983, 0, 3,
                                                                       389688, 317472, 389733,
                                                                       230964, 231048, 319308,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 392118, 0, 3,
                                                                       389733, 317508, 389778,
                                                                       231048, 231132, 319416,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 392253, 0, 3,
                                                                       389778, 317544, 389823,
                                                                       231132, 231216, 319524,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 392388, 0, 3,
                                                                       389823, 317580, 389868,
                                                                       231216, 231300, 319632,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 392523, 0, 3,
                                                                       389868, 317616, 389913,
                                                                       231300, 231384, 319740,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 392658, 0, 3,
                                                                       389958, 317688, 390093,
                                                                       231552, 231720, 319848,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 392928, 0, 3,
                                                                       390093, 317796, 390228,
                                                                       231720, 231888, 320064,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 393198, 0, 3,
                                                                       390228, 317904, 390363,
                                                                       231888, 232056, 320280,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 393468, 0, 3,
                                                                       390363, 318012, 390498,
                                                                       232056, 232224, 320496,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 393738, 0, 3,
                                                                       390498, 318120, 390633,
                                                                       232224, 232392, 320712,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 394008, 0, 3,
                                                                       390633, 318228, 390768,
                                                                       232392, 232560, 320928,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 394278, 0, 3,
                                                                       390768, 318336, 390903,
                                                                       232560, 232728, 321144,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 394548, 0, 3,
                                                                       390903, 318444, 391038,
                                                                       232728, 232896, 321360,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 394818, 0, 3,
                                                                       391038, 318552, 391173,
                                                                       232896, 233064, 321576,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 395088, 0, 3,
                                                                       391308, 318768, 391443,
                                                                       233400, 233568, 321792,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 395358, 0, 3,
                                                                       391443, 318876, 391578,
                                                                       233568, 233736, 322008,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 395628, 0, 3,
                                                                       391578, 318984, 391713,
                                                                       233736, 233904, 322224,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 395898, 0, 3,
                                                                       391713, 319092, 391848,
                                                                       233904, 234072, 322440,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 396168, 0, 3,
                                                                       391848, 319200, 391983,
                                                                       234072, 234240, 322656,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 396438, 0, 3,
                                                                       391983, 319308, 392118,
                                                                       234240, 234408, 322872,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 396708, 0, 3,
                                                                       392118, 319416, 392253,
                                                                       234408, 234576, 323088,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 396978, 0, 3,
                                                                       392253, 319524, 392388,
                                                                       234576, 234744, 323304,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 397248, 0, 3,
                                                                       392388, 319632, 392523,
                                                                       234744, 234912, 323520,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 397518, 0, 3,
                                                                       392658, 319848, 392928,
                                                                       235248, 235528, 323736,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 397968, 0, 3,
                                                                       392928, 320064, 393198,
                                                                       235528, 235808, 324096,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 398418, 0, 3,
                                                                       393198, 320280, 393468,
                                                                       235808, 236088, 324456,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 398868, 0, 3,
                                                                       393468, 320496, 393738,
                                                                       236088, 236368, 324816,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 399318, 0, 3,
                                                                       393738, 320712, 394008,
                                                                       236368, 236648, 325176,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 399768, 0, 3,
                                                                       394008, 320928, 394278,
                                                                       236648, 236928, 325536,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 400218, 0, 3,
                                                                       394278, 321144, 394548,
                                                                       236928, 237208, 325896,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 400668, 0, 3,
                                                                       394548, 321360, 394818,
                                                                       237208, 237488, 326256,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 401118, 0, 3,
                                                                       395088, 321792, 395358,
                                                                       238048, 238328, 326616,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 401568, 0, 3,
                                                                       395358, 322008, 395628,
                                                                       238328, 238608, 326976,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 402018, 0, 3,
                                                                       395628, 322224, 395898,
                                                                       238608, 238888, 327336,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 402468, 0, 3,
                                                                       395898, 322440, 396168,
                                                                       238888, 239168, 327696,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 402918, 0, 3,
                                                                       396168, 322656, 396438,
                                                                       239168, 239448, 328056,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 403368, 0, 3,
                                                                       396438, 322872, 396708,
                                                                       239448, 239728, 328416,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 403818, 0, 3,
                                                                       396708, 323088, 396978,
                                                                       239728, 240008, 328776,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 404268, 0, 3,
                                                                       396978, 323304, 397248,
                                                                       240008, 240288, 329136,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 404718, 0, 3,
                                                                       397518, 323736, 397968,
                                                                       240848, 241268, 329496,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 405393, 0, 3,
                                                                       397968, 324096, 398418,
                                                                       241268, 241688, 330036,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 406068, 0, 3,
                                                                       398418, 324456, 398868,
                                                                       241688, 242108, 330576,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 406743, 0, 3,
                                                                       398868, 324816, 399318,
                                                                       242108, 242528, 331116,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 407418, 0, 3,
                                                                       399318, 325176, 399768,
                                                                       242528, 242948, 331656,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 408093, 0, 3,
                                                                       399768, 325536, 400218,
                                                                       242948, 243368, 332196,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 408768, 0, 3,
                                                                       400218, 325896, 400668,
                                                                       243368, 243788, 332736,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 409443, 0, 3,
                                                                       401118, 326616, 401568,
                                                                       244628, 245048, 333276,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 410118, 0, 3,
                                                                       401568, 326976, 402018,
                                                                       245048, 245468, 333816,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 410793, 0, 3,
                                                                       402018, 327336, 402468,
                                                                       245468, 245888, 334356,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 411468, 0, 3,
                                                                       402468, 327696, 402918,
                                                                       245888, 246308, 334896,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 412143, 0, 3,
                                                                       402918, 328056, 403368,
                                                                       246308, 246728, 335436,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 412818, 0, 3,
                                                                       403368, 328416, 403818,
                                                                       246728, 247148, 335976,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 413493, 0, 3,
                                                                       403818, 328776, 404268,
                                                                       247148, 247568, 336516,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 414168, 0, 3,
                                                                       404718, 329496, 405393,
                                                                       248408, 248996, 337056,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 415113, 0, 3,
                                                                       405393, 330036, 406068,
                                                                       248996, 249584, 337812,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 416058, 0, 3,
                                                                       406068, 330576, 406743,
                                                                       249584, 250172, 338568,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 417003, 0, 3,
                                                                       406743, 331116, 407418,
                                                                       250172, 250760, 339324,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 417948, 0, 3,
                                                                       407418, 331656, 408093,
                                                                       250760, 251348, 340080,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 418893, 0, 3,
                                                                       408093, 332196, 408768,
                                                                       251348, 251936, 340836,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 419838, 0, 3,
                                                                       409443, 333276, 410118,
                                                                       253112, 253700, 341592,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 420783, 0, 3,
                                                                       410118, 333816, 410793,
                                                                       253700, 254288, 342348,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 421728, 0, 3,
                                                                       410793, 334356, 411468,
                                                                       254288, 254876, 343104,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 422673, 0, 3,
                                                                       411468, 334896, 412143,
                                                                       254876, 255464, 343860,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 423618, 0, 3,
                                                                       412143, 335436, 412818,
                                                                       255464, 256052, 344616,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 424563, 0, 3,
                                                                       412818, 335976, 413493,
                                                                       256052, 256640, 345372,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 425508, 0, 3,
                                                                       414168, 337056, 415113,
                                                                       257816, 258600, 346128,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 426768, 0, 3,
                                                                       415113, 337812, 416058,
                                                                       258600, 259384, 347136,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 428028, 0, 3,
                                                                       416058, 338568, 417003,
                                                                       259384, 260168, 348144,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 429288, 0, 3,
                                                                       417003, 339324, 417948,
                                                                       260168, 260952, 349152,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 430548, 0, 3,
                                                                       417948, 340080, 418893,
                                                                       260952, 261736, 350160,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 431808, 0, 3,
                                                                       419838, 341592, 420783,
                                                                       263304, 264088, 351168,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 433068, 0, 3,
                                                                       420783, 342348, 421728,
                                                                       264088, 264872, 352176,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 434328, 0, 3,
                                                                       421728, 343104, 422673,
                                                                       264872, 265656, 353184,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 435588, 0, 3,
                                                                       422673, 343860, 423618,
                                                                       265656, 266440, 354192,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 436848, 0, 3,
                                                                       423618, 344616, 424563,
                                                                       266440, 267224, 355200,
                                                                       ncols, gamma, p, q);

                    compute_prim_skl_three_center_electron_repulsion_0(buffer, 438108, 0, 3,
                                                                       425508, 346128, 426768,
                                                                       268792, 269800, 356208,
                                                                       ncols, gamma, p, q);

                    compute_prim_skl_three_center_electron_repulsion_0(buffer, 439728, 0, 3,
                                                                       426768, 347136, 428028,
                                                                       269800, 270808, 357504,
                                                                       ncols, gamma, p, q);

                    compute_prim_skl_three_center_electron_repulsion_0(buffer, 441348, 0, 3,
                                                                       428028, 348144, 429288,
                                                                       270808, 271816, 358800,
                                                                       ncols, gamma, p, q);

                    compute_prim_skl_three_center_electron_repulsion_0(buffer, 442968, 0, 3,
                                                                       429288, 349152, 430548,
                                                                       271816, 272824, 360096,
                                                                       ncols, gamma, p, q);

                    compute_prim_skl_three_center_electron_repulsion_0(buffer, 444588, 0, 3,
                                                                       431808, 351168, 433068,
                                                                       274840, 275848, 361392,
                                                                       ncols, gamma, p, q);

                    compute_prim_skl_three_center_electron_repulsion_0(buffer, 446208, 0, 3,
                                                                       433068, 352176, 434328,
                                                                       275848, 276856, 362688,
                                                                       ncols, gamma, p, q);

                    compute_prim_skl_three_center_electron_repulsion_0(buffer, 447828, 0, 3,
                                                                       434328, 353184, 435588,
                                                                       276856, 277864, 363984,
                                                                       ncols, gamma, p, q);

                    compute_prim_skl_three_center_electron_repulsion_0(buffer, 449448, 0, 3,
                                                                       435588, 354192, 436848,
                                                                       277864, 278872, 365280,
                                                                       ncols, gamma, p, q);

                    compute_prim_sll_three_center_electron_repulsion_0(buffer, 451068, 0, 3,
                                                                       438108, 356208, 439728,
                                                                       280888, 282148, 366576,
                                                                       ncols, gamma, p, q);

                    compute_prim_sll_three_center_electron_repulsion_0(buffer, 453093, 0, 3,
                                                                       439728, 357504, 441348,
                                                                       282148, 283408, 368196,
                                                                       ncols, gamma, p, q);

                    compute_prim_sll_three_center_electron_repulsion_0(buffer, 455118, 0, 3,
                                                                       441348, 358800, 442968,
                                                                       283408, 284668, 369816,
                                                                       ncols, gamma, p, q);

                    compute_prim_sll_three_center_electron_repulsion_0(buffer, 457143, 0, 3,
                                                                       444588, 361392, 446208,
                                                                       287188, 288448, 371436,
                                                                       ncols, gamma, p, q);

                    compute_prim_sll_three_center_electron_repulsion_0(buffer, 459168, 0, 3,
                                                                       446208, 362688, 447828,
                                                                       288448, 289708, 373056,
                                                                       ncols, gamma, p, q);

                    compute_prim_sll_three_center_electron_repulsion_0(buffer, 461193, 0, 3,
                                                                       447828, 363984, 449448,
                                                                       289708, 290968, 374676,
                                                                       ncols, gamma, p, q);

                    compute_prim_sml_three_center_electron_repulsion_0(buffer, 463218, 0, 3,
                                                                       451068, 366576, 453093,
                                                                       293488, 295028, 376296,
                                                                       ncols, gamma, p, q);

                    compute_prim_sml_three_center_electron_repulsion_0(buffer, 465693, 0, 3,
                                                                       453093, 368196, 455118,
                                                                       295028, 296568, 378276,
                                                                       ncols, gamma, p, q);

                    compute_prim_sml_three_center_electron_repulsion_0(buffer, 468168, 0, 3,
                                                                       457143, 371436, 459168,
                                                                       299648, 301188, 380256,
                                                                       ncols, gamma, p, q);

                    compute_prim_sml_three_center_electron_repulsion_0(buffer, 470643, 0, 3,
                                                                       459168, 373056, 461193,
                                                                       301188, 302728, 382236,
                                                                       ncols, gamma, p, q);

                    compute_prim_snl_three_center_electron_repulsion_0(buffer, 473118, 0, 3,
                                                                       463218, 376296, 465693,
                                                                       305808, 307656, 384216,
                                                                       ncols, gamma, p, q);

                    compute_prim_snl_three_center_electron_repulsion_0(buffer, 476088, 0, 3,
                                                                       468168, 380256, 470643,
                                                                       311352, 313200, 386592,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 479058, 425508, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 480794, 431808, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 482530, 438108, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 484762, 444588, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 486994, 451068, 2025, ncols);

                    simdfunc::contract_primitives(buffer, 489784, 457143, 2025, ncols);

                    simdfunc::contract_primitives(buffer, 492574, 463218, 2475, ncols);

                    simdfunc::contract_primitives(buffer, 495984, 468168, 2475, ncols);

                    simdfunc::contract_primitives(buffer, 499394, 473118, 2970, ncols);

                    simdfunc::contract_primitives(buffer, 503486, 476088, 2970, ncols);
                }
            }
        }

        simdtrf::transform_l_inner(buffer, 480318, 479058, 28, 1, nmax);

        simdtrf::transform_l_inner(buffer, 482054, 480794, 28, 1, nmax);

        simdtrf::transform_l_inner(buffer, 484150, 482530, 36, 1, nmax);

        simdtrf::transform_l_inner(buffer, 486382, 484762, 36, 1, nmax);

        simdtrf::transform_l_inner(buffer, 489019, 486994, 45, 1, nmax);

        simdtrf::transform_l_inner(buffer, 491809, 489784, 45, 1, nmax);

        simdtrf::transform_l_inner(buffer, 495049, 492574, 55, 1, nmax);

        simdtrf::transform_l_inner(buffer, 498459, 495984, 55, 1, nmax);

        simdtrf::transform_l_inner(buffer, 502364, 499394, 66, 1, nmax);

        simdtrf::transform_l_inner(buffer, 506456, 503486, 66, 1, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 507578, 480318, 484150, 17, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 509006, 482054, 486382, 17, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 510434, 484150, 489019, 17, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 512270, 486382, 491809, 17, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 514106, 489019, 495049, 17, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 516401, 491809, 498459, 17, nmax);

        simdtrf::compute_hrr_pm(buffer, coordinates, 518696, 495049, 502364, 17, nmax);

        simdtrf::compute_hrr_pm(buffer, coordinates, 521501, 498459, 506456, 17, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 524306, 507578, 510434, 17, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 527162, 509006, 512270, 17, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 530018, 510434, 514106, 17, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 533690, 512270, 516401, 17, nmax);

        simdtrf::compute_hrr_dl(buffer, coordinates, 537362, 514106, 518696, 17, nmax);

        simdtrf::compute_hrr_dl(buffer, coordinates, 541952, 516401, 521501, 17, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 546542, 524306, 530018, 17, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 551302, 527162, 533690, 17, nmax);

        simdtrf::compute_hrr_fk(buffer, coordinates, 556062, 530018, 537362, 17, nmax);

        simdtrf::compute_hrr_fk(buffer, coordinates, 562182, 533690, 541952, 17, nmax);

        simdtrf::compute_hrr_gi(buffer, coordinates, 568302, 546542, 556062, 17, nmax);

        simdtrf::compute_hrr_gi(buffer, coordinates, 575442, 551302, 562182, 17, nmax);

        simdtrf::transform_i_inner(buffer, 582582, 575442, 15, 17, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 582582, 221, nmax);

        simdtrf::transform_i_inner(buffer, 582582, 568302, 15, 17, nmax);

        simdtrf::transform_g_outer(values + 1989 * nvalues + n * npairs, nvalues, buffer, 582582,
                                   221, nmax);
    }

    for (size_t m = 0; m < 3978; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
