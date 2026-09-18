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


#include "SimdThreeCenterElectronRepulsionRsRecHIL.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSOD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOS.hpp"
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
#include "SimdTransformL.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_hil_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_hil_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 829179, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 4862 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 829179, 639738, 36866, dimensions);

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

                    simdfunc::compute_full_t3c_erf_boys_function(buffer, coordinates, 6, 3, 19,
                                                                 ncols, fj, i * nprim_b + j, fq,
                                                                 omega);

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 27, 3, 19,
                                                             ncols, fj, i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 48, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 51, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 54, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 57, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 60, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 63, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 66, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 69, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 72, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 75, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 78, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 81, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 84, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 87, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 90, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 93, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 96, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 99, 0, 3, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 102, 0, 3, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 105, 0, 3, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 108, 0, 3, 29, 30,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 111, 0, 3, 30, 31,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 114, 0, 3, 31, 32,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 117, 0, 3, 32, 33,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 120, 0, 3, 33, 34,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 123, 0, 3, 34, 35,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 126, 0, 3, 35, 36,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 129, 0, 3, 36, 37,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 132, 0, 3, 37, 38,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 135, 0, 3, 38, 39,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 138, 0, 3, 39, 40,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 141, 0, 3, 40, 41,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 144, 0, 3, 41, 42,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 147, 0, 3, 42, 43,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 150, 0, 3, 43, 44,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 153, 0, 3, 44, 45,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 156, 0, 3, 45, 46,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 159, 0, 3, 46, 47,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 162, 0, 3, 7, 8,
                                                                       48, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 168, 0, 3, 8, 9,
                                                                       51, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 174, 0, 3, 9, 10,
                                                                       54, 57, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 180, 0, 3, 10, 11,
                                                                       57, 60, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 186, 0, 3, 11, 12,
                                                                       60, 63, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 192, 0, 3, 12, 13,
                                                                       63, 66, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 198, 0, 3, 13, 14,
                                                                       66, 69, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 204, 0, 3, 14, 15,
                                                                       69, 72, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 210, 0, 3, 15, 16,
                                                                       72, 75, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 216, 0, 3, 16, 17,
                                                                       75, 78, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 222, 0, 3, 17, 18,
                                                                       78, 81, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 228, 0, 3, 18, 19,
                                                                       81, 84, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 234, 0, 3, 19, 20,
                                                                       84, 87, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 240, 0, 3, 20, 21,
                                                                       87, 90, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 246, 0, 3, 21, 22,
                                                                       90, 93, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 252, 0, 3, 22, 23,
                                                                       93, 96, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 258, 0, 3, 23, 24,
                                                                       96, 99, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 264, 0, 3, 24, 25,
                                                                       99, 102, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 270, 0, 3, 28, 29,
                                                                       105, 108, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 276, 0, 3, 29, 30,
                                                                       108, 111, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 282, 0, 3, 30, 31,
                                                                       111, 114, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 288, 0, 3, 31, 32,
                                                                       114, 117, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 294, 0, 3, 32, 33,
                                                                       117, 120, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 300, 0, 3, 33, 34,
                                                                       120, 123, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 306, 0, 3, 34, 35,
                                                                       123, 126, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 312, 0, 3, 35, 36,
                                                                       126, 129, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 318, 0, 3, 36, 37,
                                                                       129, 132, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 324, 0, 3, 37, 38,
                                                                       132, 135, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 330, 0, 3, 38, 39,
                                                                       135, 138, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 336, 0, 3, 39, 40,
                                                                       138, 141, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 342, 0, 3, 40, 41,
                                                                       141, 144, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 348, 0, 3, 41, 42,
                                                                       144, 147, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 354, 0, 3, 42, 43,
                                                                       147, 150, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 360, 0, 3, 43, 44,
                                                                       150, 153, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 366, 0, 3, 44, 45,
                                                                       153, 156, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 372, 0, 3, 45, 46,
                                                                       156, 159, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 378, 0, 3, 48, 51,
                                                                       162, 168, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 388, 0, 3, 51, 54,
                                                                       168, 174, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 398, 0, 3, 54, 57,
                                                                       174, 180, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 408, 0, 3, 57, 60,
                                                                       180, 186, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 418, 0, 3, 60, 63,
                                                                       186, 192, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 428, 0, 3, 63, 66,
                                                                       192, 198, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 438, 0, 3, 66, 69,
                                                                       198, 204, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 448, 0, 3, 69, 72,
                                                                       204, 210, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 458, 0, 3, 72, 75,
                                                                       210, 216, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 468, 0, 3, 75, 78,
                                                                       216, 222, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 478, 0, 3, 78, 81,
                                                                       222, 228, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 488, 0, 3, 81, 84,
                                                                       228, 234, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 498, 0, 3, 84, 87,
                                                                       234, 240, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 508, 0, 3, 87, 90,
                                                                       240, 246, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 518, 0, 3, 90, 93,
                                                                       246, 252, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 528, 0, 3, 93, 96,
                                                                       252, 258, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 538, 0, 3, 96, 99,
                                                                       258, 264, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 548, 0, 3, 105,
                                                                       108, 270, 276, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 558, 0, 3, 108,
                                                                       111, 276, 282, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 568, 0, 3, 111,
                                                                       114, 282, 288, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 578, 0, 3, 114,
                                                                       117, 288, 294, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 588, 0, 3, 117,
                                                                       120, 294, 300, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 598, 0, 3, 120,
                                                                       123, 300, 306, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 608, 0, 3, 123,
                                                                       126, 306, 312, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 618, 0, 3, 126,
                                                                       129, 312, 318, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 628, 0, 3, 129,
                                                                       132, 318, 324, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 638, 0, 3, 132,
                                                                       135, 324, 330, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 648, 0, 3, 135,
                                                                       138, 330, 336, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 658, 0, 3, 138,
                                                                       141, 336, 342, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 668, 0, 3, 141,
                                                                       144, 342, 348, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 678, 0, 3, 144,
                                                                       147, 348, 354, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 688, 0, 3, 147,
                                                                       150, 354, 360, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 698, 0, 3, 150,
                                                                       153, 360, 366, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 708, 0, 3, 153,
                                                                       156, 366, 372, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 718, 0, 3, 162,
                                                                       168, 378, 388, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 733, 0, 3, 168,
                                                                       174, 388, 398, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 748, 0, 3, 174,
                                                                       180, 398, 408, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 763, 0, 3, 180,
                                                                       186, 408, 418, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 778, 0, 3, 186,
                                                                       192, 418, 428, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 793, 0, 3, 192,
                                                                       198, 428, 438, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 808, 0, 3, 198,
                                                                       204, 438, 448, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 823, 0, 3, 204,
                                                                       210, 448, 458, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 838, 0, 3, 210,
                                                                       216, 458, 468, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 853, 0, 3, 216,
                                                                       222, 468, 478, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 868, 0, 3, 222,
                                                                       228, 478, 488, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 883, 0, 3, 228,
                                                                       234, 488, 498, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 898, 0, 3, 234,
                                                                       240, 498, 508, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 913, 0, 3, 240,
                                                                       246, 508, 518, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 928, 0, 3, 246,
                                                                       252, 518, 528, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 943, 0, 3, 252,
                                                                       258, 528, 538, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 958, 0, 3, 270,
                                                                       276, 548, 558, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 973, 0, 3, 276,
                                                                       282, 558, 568, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 988, 0, 3, 282,
                                                                       288, 568, 578, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 1003, 0, 3, 288,
                                                                       294, 578, 588, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 1018, 0, 3, 294,
                                                                       300, 588, 598, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 1033, 0, 3, 300,
                                                                       306, 598, 608, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 1048, 0, 3, 306,
                                                                       312, 608, 618, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 1063, 0, 3, 312,
                                                                       318, 618, 628, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 1078, 0, 3, 318,
                                                                       324, 628, 638, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 1093, 0, 3, 324,
                                                                       330, 638, 648, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 1108, 0, 3, 330,
                                                                       336, 648, 658, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 1123, 0, 3, 336,
                                                                       342, 658, 668, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 1138, 0, 3, 342,
                                                                       348, 668, 678, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 1153, 0, 3, 348,
                                                                       354, 678, 688, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 1168, 0, 3, 354,
                                                                       360, 688, 698, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 1183, 0, 3, 360,
                                                                       366, 698, 708, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1198, 0, 3, 378,
                                                                       388, 718, 733, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1219, 0, 3, 388,
                                                                       398, 733, 748, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1240, 0, 3, 398,
                                                                       408, 748, 763, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1261, 0, 3, 408,
                                                                       418, 763, 778, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1282, 0, 3, 418,
                                                                       428, 778, 793, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1303, 0, 3, 428,
                                                                       438, 793, 808, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1324, 0, 3, 438,
                                                                       448, 808, 823, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1345, 0, 3, 448,
                                                                       458, 823, 838, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1366, 0, 3, 458,
                                                                       468, 838, 853, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1387, 0, 3, 468,
                                                                       478, 853, 868, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1408, 0, 3, 478,
                                                                       488, 868, 883, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1429, 0, 3, 488,
                                                                       498, 883, 898, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1450, 0, 3, 498,
                                                                       508, 898, 913, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1471, 0, 3, 508,
                                                                       518, 913, 928, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1492, 0, 3, 518,
                                                                       528, 928, 943, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1513, 0, 3, 548,
                                                                       558, 958, 973, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1534, 0, 3, 558,
                                                                       568, 973, 988, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1555, 0, 3, 568,
                                                                       578, 988, 1003, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1576, 0, 3, 578,
                                                                       588, 1003, 1018, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1597, 0, 3, 588,
                                                                       598, 1018, 1033, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1618, 0, 3, 598,
                                                                       608, 1033, 1048, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1639, 0, 3, 608,
                                                                       618, 1048, 1063, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1660, 0, 3, 618,
                                                                       628, 1063, 1078, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1681, 0, 3, 628,
                                                                       638, 1078, 1093, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1702, 0, 3, 638,
                                                                       648, 1093, 1108, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1723, 0, 3, 648,
                                                                       658, 1108, 1123, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1744, 0, 3, 658,
                                                                       668, 1123, 1138, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1765, 0, 3, 668,
                                                                       678, 1138, 1153, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1786, 0, 3, 678,
                                                                       688, 1153, 1168, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1807, 0, 3, 688,
                                                                       698, 1168, 1183, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1828, 0, 3, 718,
                                                                       733, 1198, 1219, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1856, 0, 3, 733,
                                                                       748, 1219, 1240, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1884, 0, 3, 748,
                                                                       763, 1240, 1261, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1912, 0, 3, 763,
                                                                       778, 1261, 1282, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1940, 0, 3, 778,
                                                                       793, 1282, 1303, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1968, 0, 3, 793,
                                                                       808, 1303, 1324, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1996, 0, 3, 808,
                                                                       823, 1324, 1345, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2024, 0, 3, 823,
                                                                       838, 1345, 1366, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2052, 0, 3, 838,
                                                                       853, 1366, 1387, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2080, 0, 3, 853,
                                                                       868, 1387, 1408, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2108, 0, 3, 868,
                                                                       883, 1408, 1429, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2136, 0, 3, 883,
                                                                       898, 1429, 1450, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2164, 0, 3, 898,
                                                                       913, 1450, 1471, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2192, 0, 3, 913,
                                                                       928, 1471, 1492, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2220, 0, 3, 958,
                                                                       973, 1513, 1534, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2248, 0, 3, 973,
                                                                       988, 1534, 1555, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2276, 0, 3, 988,
                                                                       1003, 1555, 1576, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2304, 0, 3, 1003,
                                                                       1018, 1576, 1597, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2332, 0, 3, 1018,
                                                                       1033, 1597, 1618, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2360, 0, 3, 1033,
                                                                       1048, 1618, 1639, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2388, 0, 3, 1048,
                                                                       1063, 1639, 1660, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2416, 0, 3, 1063,
                                                                       1078, 1660, 1681, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2444, 0, 3, 1078,
                                                                       1093, 1681, 1702, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2472, 0, 3, 1093,
                                                                       1108, 1702, 1723, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2500, 0, 3, 1108,
                                                                       1123, 1723, 1744, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2528, 0, 3, 1123,
                                                                       1138, 1744, 1765, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2556, 0, 3, 1138,
                                                                       1153, 1765, 1786, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2584, 0, 3, 1153,
                                                                       1168, 1786, 1807, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2612, 0, 3, 1198,
                                                                       1219, 1828, 1856, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2648, 0, 3, 1219,
                                                                       1240, 1856, 1884, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2684, 0, 3, 1240,
                                                                       1261, 1884, 1912, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2720, 0, 3, 1261,
                                                                       1282, 1912, 1940, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2756, 0, 3, 1282,
                                                                       1303, 1940, 1968, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2792, 0, 3, 1303,
                                                                       1324, 1968, 1996, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2828, 0, 3, 1324,
                                                                       1345, 1996, 2024, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2864, 0, 3, 1345,
                                                                       1366, 2024, 2052, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2900, 0, 3, 1366,
                                                                       1387, 2052, 2080, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2936, 0, 3, 1387,
                                                                       1408, 2080, 2108, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2972, 0, 3, 1408,
                                                                       1429, 2108, 2136, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 3008, 0, 3, 1429,
                                                                       1450, 2136, 2164, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 3044, 0, 3, 1450,
                                                                       1471, 2164, 2192, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 3080, 0, 3, 1513,
                                                                       1534, 2220, 2248, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 3116, 0, 3, 1534,
                                                                       1555, 2248, 2276, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 3152, 0, 3, 1555,
                                                                       1576, 2276, 2304, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 3188, 0, 3, 1576,
                                                                       1597, 2304, 2332, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 3224, 0, 3, 1597,
                                                                       1618, 2332, 2360, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 3260, 0, 3, 1618,
                                                                       1639, 2360, 2388, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 3296, 0, 3, 1639,
                                                                       1660, 2388, 2416, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 3332, 0, 3, 1660,
                                                                       1681, 2416, 2444, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 3368, 0, 3, 1681,
                                                                       1702, 2444, 2472, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 3404, 0, 3, 1702,
                                                                       1723, 2472, 2500, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 3440, 0, 3, 1723,
                                                                       1744, 2500, 2528, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 3476, 0, 3, 1744,
                                                                       1765, 2528, 2556, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 3512, 0, 3, 1765,
                                                                       1786, 2556, 2584, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3548, 0, 3, 1828,
                                                                       1856, 2612, 2648, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3593, 0, 3, 1856,
                                                                       1884, 2648, 2684, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3638, 0, 3, 1884,
                                                                       1912, 2684, 2720, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3683, 0, 3, 1912,
                                                                       1940, 2720, 2756, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3728, 0, 3, 1940,
                                                                       1968, 2756, 2792, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3773, 0, 3, 1968,
                                                                       1996, 2792, 2828, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3818, 0, 3, 1996,
                                                                       2024, 2828, 2864, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3863, 0, 3, 2024,
                                                                       2052, 2864, 2900, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3908, 0, 3, 2052,
                                                                       2080, 2900, 2936, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3953, 0, 3, 2080,
                                                                       2108, 2936, 2972, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3998, 0, 3, 2108,
                                                                       2136, 2972, 3008, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 4043, 0, 3, 2136,
                                                                       2164, 3008, 3044, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 4088, 0, 3, 2220,
                                                                       2248, 3080, 3116, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 4133, 0, 3, 2248,
                                                                       2276, 3116, 3152, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 4178, 0, 3, 2276,
                                                                       2304, 3152, 3188, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 4223, 0, 3, 2304,
                                                                       2332, 3188, 3224, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 4268, 0, 3, 2332,
                                                                       2360, 3224, 3260, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 4313, 0, 3, 2360,
                                                                       2388, 3260, 3296, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 4358, 0, 3, 2388,
                                                                       2416, 3296, 3332, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 4403, 0, 3, 2416,
                                                                       2444, 3332, 3368, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 4448, 0, 3, 2444,
                                                                       2472, 3368, 3404, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 4493, 0, 3, 2472,
                                                                       2500, 3404, 3440, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 4538, 0, 3, 2500,
                                                                       2528, 3440, 3476, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 4583, 0, 3, 2528,
                                                                       2556, 3476, 3512, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4628, 0, 3, 2612,
                                                                       2648, 3548, 3593, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4683, 0, 3, 2648,
                                                                       2684, 3593, 3638, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4738, 0, 3, 2684,
                                                                       2720, 3638, 3683, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4793, 0, 3, 2720,
                                                                       2756, 3683, 3728, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4848, 0, 3, 2756,
                                                                       2792, 3728, 3773, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4903, 0, 3, 2792,
                                                                       2828, 3773, 3818, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4958, 0, 3, 2828,
                                                                       2864, 3818, 3863, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 5013, 0, 3, 2864,
                                                                       2900, 3863, 3908, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 5068, 0, 3, 2900,
                                                                       2936, 3908, 3953, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 5123, 0, 3, 2936,
                                                                       2972, 3953, 3998, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 5178, 0, 3, 2972,
                                                                       3008, 3998, 4043, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 5233, 0, 3, 3080,
                                                                       3116, 4088, 4133, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 5288, 0, 3, 3116,
                                                                       3152, 4133, 4178, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 5343, 0, 3, 3152,
                                                                       3188, 4178, 4223, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 5398, 0, 3, 3188,
                                                                       3224, 4223, 4268, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 5453, 0, 3, 3224,
                                                                       3260, 4268, 4313, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 5508, 0, 3, 3260,
                                                                       3296, 4313, 4358, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 5563, 0, 3, 3296,
                                                                       3332, 4358, 4403, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 5618, 0, 3, 3332,
                                                                       3368, 4403, 4448, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 5673, 0, 3, 3368,
                                                                       3404, 4448, 4493, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 5728, 0, 3, 3404,
                                                                       3440, 4493, 4538, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 5783, 0, 3, 3440,
                                                                       3476, 4538, 4583, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 5838, 0, 3, 3548,
                                                                       3593, 4628, 4683, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 5904, 0, 3, 3593,
                                                                       3638, 4683, 4738, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 5970, 0, 3, 3638,
                                                                       3683, 4738, 4793, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 6036, 0, 3, 3683,
                                                                       3728, 4793, 4848, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 6102, 0, 3, 3728,
                                                                       3773, 4848, 4903, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 6168, 0, 3, 3773,
                                                                       3818, 4903, 4958, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 6234, 0, 3, 3818,
                                                                       3863, 4958, 5013, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 6300, 0, 3, 3863,
                                                                       3908, 5013, 5068, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 6366, 0, 3, 3908,
                                                                       3953, 5068, 5123, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 6432, 0, 3, 3953,
                                                                       3998, 5123, 5178, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 6498, 0, 3, 4088,
                                                                       4133, 5233, 5288, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 6564, 0, 3, 4133,
                                                                       4178, 5288, 5343, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 6630, 0, 3, 4178,
                                                                       4223, 5343, 5398, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 6696, 0, 3, 4223,
                                                                       4268, 5398, 5453, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 6762, 0, 3, 4268,
                                                                       4313, 5453, 5508, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 6828, 0, 3, 4313,
                                                                       4358, 5508, 5563, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 6894, 0, 3, 4358,
                                                                       4403, 5563, 5618, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 6960, 0, 3, 4403,
                                                                       4448, 5618, 5673, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 7026, 0, 3, 4448,
                                                                       4493, 5673, 5728, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 7092, 0, 3, 4493,
                                                                       4538, 5728, 5783, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 7158, 0, 3, 4628,
                                                                       4683, 5838, 5904, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 7236, 0, 3, 4683,
                                                                       4738, 5904, 5970, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 7314, 0, 3, 4738,
                                                                       4793, 5970, 6036, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 7392, 0, 3, 4793,
                                                                       4848, 6036, 6102, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 7470, 0, 3, 4848,
                                                                       4903, 6102, 6168, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 7548, 0, 3, 4903,
                                                                       4958, 6168, 6234, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 7626, 0, 3, 4958,
                                                                       5013, 6234, 6300, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 7704, 0, 3, 5013,
                                                                       5068, 6300, 6366, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 7782, 0, 3, 5068,
                                                                       5123, 6366, 6432, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 7860, 0, 3, 5233,
                                                                       5288, 6498, 6564, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 7938, 0, 3, 5288,
                                                                       5343, 6564, 6630, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 8016, 0, 3, 5343,
                                                                       5398, 6630, 6696, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 8094, 0, 3, 5398,
                                                                       5453, 6696, 6762, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 8172, 0, 3, 5453,
                                                                       5508, 6762, 6828, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 8250, 0, 3, 5508,
                                                                       5563, 6828, 6894, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 8328, 0, 3, 5563,
                                                                       5618, 6894, 6960, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 8406, 0, 3, 5618,
                                                                       5673, 6960, 7026, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 8484, 0, 3, 5673,
                                                                       5728, 7026, 7092, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 8562, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 8565, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 8568, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 8571, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 8574, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 8577, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 8580, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 8583, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 8586, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 8589, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 8592, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 8595, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 8598, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 8601, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 8604, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 8607, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 8610, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 8613, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 8616, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 8619, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 8622, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 8625, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 8628, 3, 34,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 8631, 3, 35,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 8634, 3, 36,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 8637, 3, 37,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 8640, 3, 38,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 8643, 3, 39,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 8646, 3, 40,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 8649, 3, 41,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 8652, 3, 42,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 8655, 3, 43,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 8658, 3, 44,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 8661, 3, 45,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 8664, 3, 46,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 8667, 3, 47,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 8670, 3, 9, 54,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 8679, 3, 10, 57,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 8688, 3, 11, 60,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 8697, 3, 12, 63,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 8706, 3, 13, 66,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 8715, 3, 14, 69,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 8724, 3, 15, 72,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 8733, 3, 16, 75,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 8742, 3, 17, 78,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 8751, 3, 18, 81,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 8760, 3, 19, 84,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 8769, 3, 20, 87,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 8778, 3, 21, 90,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 8787, 3, 22, 93,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 8796, 3, 23, 96,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 8805, 3, 24, 99,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 8814, 3, 25, 102,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 8823, 3, 30, 111,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 8832, 3, 31, 114,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 8841, 3, 32, 117,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 8850, 3, 33, 120,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 8859, 3, 34, 123,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 8868, 3, 35, 126,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 8877, 3, 36, 129,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 8886, 3, 37, 132,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 8895, 3, 38, 135,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 8904, 3, 39, 138,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 8913, 3, 40, 141,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 8922, 3, 41, 144,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 8931, 3, 42, 147,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 8940, 3, 43, 150,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 8949, 3, 44, 153,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 8958, 3, 45, 156,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 8967, 3, 46, 159,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 8976, 3, 54, 174,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 8994, 3, 57, 180,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 9012, 3, 60, 186,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 9030, 3, 63, 192,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 9048, 3, 66, 198,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 9066, 3, 69, 204,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 9084, 3, 72, 210,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 9102, 3, 75, 216,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 9120, 3, 78, 222,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 9138, 3, 81, 228,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 9156, 3, 84, 234,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 9174, 3, 87, 240,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 9192, 3, 90, 246,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 9210, 3, 93, 252,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 9228, 3, 96, 258,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 9246, 3, 99, 264,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 9264, 3, 111, 282,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 9282, 3, 114, 288,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 9300, 3, 117, 294,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 9318, 3, 120, 300,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 9336, 3, 123, 306,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 9354, 3, 126, 312,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 9372, 3, 129, 318,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 9390, 3, 132, 324,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 9408, 3, 135, 330,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 9426, 3, 138, 336,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 9444, 3, 141, 342,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 9462, 3, 144, 348,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 9480, 3, 147, 354,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 9498, 3, 150, 360,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 9516, 3, 153, 366,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 9534, 3, 156, 372,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 9552, 3, 174, 398,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 9582, 3, 180, 408,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 9612, 3, 186, 418,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 9642, 3, 192, 428,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 9672, 3, 198, 438,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 9702, 3, 204, 448,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 9732, 3, 210, 458,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 9762, 3, 216, 468,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 9792, 3, 222, 478,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 9822, 3, 228, 488,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 9852, 3, 234, 498,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 9882, 3, 240, 508,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 9912, 3, 246, 518,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 9942, 3, 252, 528,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 9972, 3, 258, 538,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 10002, 3, 282,
                                                                       568, ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 10032, 3, 288,
                                                                       578, ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 10062, 3, 294,
                                                                       588, ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 10092, 3, 300,
                                                                       598, ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 10122, 3, 306,
                                                                       608, ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 10152, 3, 312,
                                                                       618, ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 10182, 3, 318,
                                                                       628, ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 10212, 3, 324,
                                                                       638, ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 10242, 3, 330,
                                                                       648, ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 10272, 3, 336,
                                                                       658, ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 10302, 3, 342,
                                                                       668, ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 10332, 3, 348,
                                                                       678, ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 10362, 3, 354,
                                                                       688, ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 10392, 3, 360,
                                                                       698, ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 10422, 3, 366,
                                                                       708, ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 10452, 3, 398,
                                                                       748, ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 10497, 3, 408,
                                                                       763, ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 10542, 3, 418,
                                                                       778, ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 10587, 3, 428,
                                                                       793, ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 10632, 3, 438,
                                                                       808, ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 10677, 3, 448,
                                                                       823, ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 10722, 3, 458,
                                                                       838, ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 10767, 3, 468,
                                                                       853, ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 10812, 3, 478,
                                                                       868, ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 10857, 3, 488,
                                                                       883, ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 10902, 3, 498,
                                                                       898, ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 10947, 3, 508,
                                                                       913, ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 10992, 3, 518,
                                                                       928, ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 11037, 3, 528,
                                                                       943, ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 11082, 3, 568,
                                                                       988, ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 11127, 3, 578,
                                                                       1003, ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 11172, 3, 588,
                                                                       1018, ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 11217, 3, 598,
                                                                       1033, ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 11262, 3, 608,
                                                                       1048, ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 11307, 3, 618,
                                                                       1063, ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 11352, 3, 628,
                                                                       1078, ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 11397, 3, 638,
                                                                       1093, ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 11442, 3, 648,
                                                                       1108, ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 11487, 3, 658,
                                                                       1123, ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 11532, 3, 668,
                                                                       1138, ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 11577, 3, 678,
                                                                       1153, ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 11622, 3, 688,
                                                                       1168, ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 11667, 3, 698,
                                                                       1183, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 11712, 3, 748,
                                                                       1240, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 11775, 3, 763,
                                                                       1261, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 11838, 3, 778,
                                                                       1282, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 11901, 3, 793,
                                                                       1303, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 11964, 3, 808,
                                                                       1324, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 12027, 3, 823,
                                                                       1345, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 12090, 3, 838,
                                                                       1366, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 12153, 3, 853,
                                                                       1387, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 12216, 3, 868,
                                                                       1408, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 12279, 3, 883,
                                                                       1429, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 12342, 3, 898,
                                                                       1450, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 12405, 3, 913,
                                                                       1471, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 12468, 3, 928,
                                                                       1492, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 12531, 3, 988,
                                                                       1555, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 12594, 3, 1003,
                                                                       1576, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 12657, 3, 1018,
                                                                       1597, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 12720, 3, 1033,
                                                                       1618, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 12783, 3, 1048,
                                                                       1639, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 12846, 3, 1063,
                                                                       1660, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 12909, 3, 1078,
                                                                       1681, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 12972, 3, 1093,
                                                                       1702, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 13035, 3, 1108,
                                                                       1723, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 13098, 3, 1123,
                                                                       1744, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 13161, 3, 1138,
                                                                       1765, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 13224, 3, 1153,
                                                                       1786, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 13287, 3, 1168,
                                                                       1807, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 13350, 3, 1240,
                                                                       1884, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 13434, 3, 1261,
                                                                       1912, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 13518, 3, 1282,
                                                                       1940, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 13602, 3, 1303,
                                                                       1968, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 13686, 3, 1324,
                                                                       1996, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 13770, 3, 1345,
                                                                       2024, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 13854, 3, 1366,
                                                                       2052, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 13938, 3, 1387,
                                                                       2080, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 14022, 3, 1408,
                                                                       2108, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 14106, 3, 1429,
                                                                       2136, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 14190, 3, 1450,
                                                                       2164, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 14274, 3, 1471,
                                                                       2192, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 14358, 3, 1555,
                                                                       2276, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 14442, 3, 1576,
                                                                       2304, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 14526, 3, 1597,
                                                                       2332, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 14610, 3, 1618,
                                                                       2360, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 14694, 3, 1639,
                                                                       2388, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 14778, 3, 1660,
                                                                       2416, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 14862, 3, 1681,
                                                                       2444, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 14946, 3, 1702,
                                                                       2472, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 15030, 3, 1723,
                                                                       2500, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 15114, 3, 1744,
                                                                       2528, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 15198, 3, 1765,
                                                                       2556, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 15282, 3, 1786,
                                                                       2584, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 15366, 3, 1884,
                                                                       2684, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 15474, 3, 1912,
                                                                       2720, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 15582, 3, 1940,
                                                                       2756, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 15690, 3, 1968,
                                                                       2792, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 15798, 3, 1996,
                                                                       2828, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 15906, 3, 2024,
                                                                       2864, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 16014, 3, 2052,
                                                                       2900, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 16122, 3, 2080,
                                                                       2936, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 16230, 3, 2108,
                                                                       2972, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 16338, 3, 2136,
                                                                       3008, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 16446, 3, 2164,
                                                                       3044, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 16554, 3, 2276,
                                                                       3152, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 16662, 3, 2304,
                                                                       3188, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 16770, 3, 2332,
                                                                       3224, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 16878, 3, 2360,
                                                                       3260, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 16986, 3, 2388,
                                                                       3296, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 17094, 3, 2416,
                                                                       3332, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 17202, 3, 2444,
                                                                       3368, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 17310, 3, 2472,
                                                                       3404, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 17418, 3, 2500,
                                                                       3440, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 17526, 3, 2528,
                                                                       3476, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 17634, 3, 2556,
                                                                       3512, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 17742, 3, 2684,
                                                                       3638, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 17877, 3, 2720,
                                                                       3683, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 18012, 3, 2756,
                                                                       3728, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 18147, 3, 2792,
                                                                       3773, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 18282, 3, 2828,
                                                                       3818, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 18417, 3, 2864,
                                                                       3863, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 18552, 3, 2900,
                                                                       3908, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 18687, 3, 2936,
                                                                       3953, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 18822, 3, 2972,
                                                                       3998, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 18957, 3, 3008,
                                                                       4043, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 19092, 3, 3152,
                                                                       4178, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 19227, 3, 3188,
                                                                       4223, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 19362, 3, 3224,
                                                                       4268, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 19497, 3, 3260,
                                                                       4313, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 19632, 3, 3296,
                                                                       4358, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 19767, 3, 3332,
                                                                       4403, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 19902, 3, 3368,
                                                                       4448, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 20037, 3, 3404,
                                                                       4493, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 20172, 3, 3440,
                                                                       4538, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 20307, 3, 3476,
                                                                       4583, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 20442, 3, 3638,
                                                                       4738, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 20607, 3, 3683,
                                                                       4793, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 20772, 3, 3728,
                                                                       4848, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 20937, 3, 3773,
                                                                       4903, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 21102, 3, 3818,
                                                                       4958, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 21267, 3, 3863,
                                                                       5013, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 21432, 3, 3908,
                                                                       5068, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 21597, 3, 3953,
                                                                       5123, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 21762, 3, 3998,
                                                                       5178, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 21927, 3, 4178,
                                                                       5343, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 22092, 3, 4223,
                                                                       5398, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 22257, 3, 4268,
                                                                       5453, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 22422, 3, 4313,
                                                                       5508, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 22587, 3, 4358,
                                                                       5563, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 22752, 3, 4403,
                                                                       5618, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 22917, 3, 4448,
                                                                       5673, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 23082, 3, 4493,
                                                                       5728, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 23247, 3, 4538,
                                                                       5783, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 23412, 3, 4738,
                                                                       5970, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 23610, 3, 4793,
                                                                       6036, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 23808, 3, 4848,
                                                                       6102, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 24006, 3, 4903,
                                                                       6168, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 24204, 3, 4958,
                                                                       6234, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 24402, 3, 5013,
                                                                       6300, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 24600, 3, 5068,
                                                                       6366, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 24798, 3, 5123,
                                                                       6432, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 24996, 3, 5343,
                                                                       6630, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 25194, 3, 5398,
                                                                       6696, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 25392, 3, 5453,
                                                                       6762, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 25590, 3, 5508,
                                                                       6828, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 25788, 3, 5563,
                                                                       6894, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 25986, 3, 5618,
                                                                       6960, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 26184, 3, 5673,
                                                                       7026, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 26382, 3, 5728,
                                                                       7092, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 26580, 3, 5970,
                                                                       7314, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 26814, 3, 6036,
                                                                       7392, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 27048, 3, 6102,
                                                                       7470, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 27282, 3, 6168,
                                                                       7548, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 27516, 3, 6234,
                                                                       7626, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 27750, 3, 6300,
                                                                       7704, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 27984, 3, 6366,
                                                                       7782, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 28218, 3, 6630,
                                                                       8016, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 28452, 3, 6696,
                                                                       8094, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 28686, 3, 6762,
                                                                       8172, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 28920, 3, 6828,
                                                                       8250, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 29154, 3, 6894,
                                                                       8328, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 29388, 3, 6960,
                                                                       8406, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 29622, 3, 7026,
                                                                       8484, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29856, 3, 7, 8,
                                                                       8562, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29862, 3, 8, 9,
                                                                       8565, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29868, 3, 9, 10,
                                                                       8568, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29874, 3, 10, 11,
                                                                       8571, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29880, 3, 11, 12,
                                                                       8574, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29886, 3, 12, 13,
                                                                       8577, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29892, 3, 13, 14,
                                                                       8580, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29898, 3, 14, 15,
                                                                       8583, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29904, 3, 15, 16,
                                                                       8586, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29910, 3, 16, 17,
                                                                       8589, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29916, 3, 17, 18,
                                                                       8592, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29922, 3, 18, 19,
                                                                       8595, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29928, 3, 19, 20,
                                                                       8598, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29934, 3, 20, 21,
                                                                       8601, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29940, 3, 21, 22,
                                                                       8604, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29946, 3, 22, 23,
                                                                       8607, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29952, 3, 23, 24,
                                                                       8610, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29958, 3, 24, 25,
                                                                       8613, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29964, 3, 28, 29,
                                                                       8616, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29970, 3, 29, 30,
                                                                       8619, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29976, 3, 30, 31,
                                                                       8622, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29982, 3, 31, 32,
                                                                       8625, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29988, 3, 32, 33,
                                                                       8628, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29994, 3, 33, 34,
                                                                       8631, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 30000, 3, 34, 35,
                                                                       8634, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 30006, 3, 35, 36,
                                                                       8637, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 30012, 3, 36, 37,
                                                                       8640, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 30018, 3, 37, 38,
                                                                       8643, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 30024, 3, 38, 39,
                                                                       8646, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 30030, 3, 39, 40,
                                                                       8649, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 30036, 3, 40, 41,
                                                                       8652, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 30042, 3, 41, 42,
                                                                       8655, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 30048, 3, 42, 43,
                                                                       8658, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 30054, 3, 43, 44,
                                                                       8661, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 30060, 3, 44, 45,
                                                                       8664, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 30066, 3, 45, 46,
                                                                       8667, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 30072, 0, 3,
                                                                       29856, 8562, 29862, 8670,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 30090, 0, 3,
                                                                       29862, 8565, 29868, 8679,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 30108, 0, 3,
                                                                       29868, 8568, 29874, 8688,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 30126, 0, 3,
                                                                       29874, 8571, 29880, 8697,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 30144, 0, 3,
                                                                       29880, 8574, 29886, 8706,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 30162, 0, 3,
                                                                       29886, 8577, 29892, 8715,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 30180, 0, 3,
                                                                       29892, 8580, 29898, 8724,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 30198, 0, 3,
                                                                       29898, 8583, 29904, 8733,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 30216, 0, 3,
                                                                       29904, 8586, 29910, 8742,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 30234, 0, 3,
                                                                       29910, 8589, 29916, 8751,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 30252, 0, 3,
                                                                       29916, 8592, 29922, 8760,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 30270, 0, 3,
                                                                       29922, 8595, 29928, 8769,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 30288, 0, 3,
                                                                       29928, 8598, 29934, 8778,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 30306, 0, 3,
                                                                       29934, 8601, 29940, 8787,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 30324, 0, 3,
                                                                       29940, 8604, 29946, 8796,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 30342, 0, 3,
                                                                       29946, 8607, 29952, 8805,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 30360, 0, 3,
                                                                       29952, 8610, 29958, 8814,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 30378, 0, 3,
                                                                       29964, 8616, 29970, 8823,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 30396, 0, 3,
                                                                       29970, 8619, 29976, 8832,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 30414, 0, 3,
                                                                       29976, 8622, 29982, 8841,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 30432, 0, 3,
                                                                       29982, 8625, 29988, 8850,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 30450, 0, 3,
                                                                       29988, 8628, 29994, 8859,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 30468, 0, 3,
                                                                       29994, 8631, 30000, 8868,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 30486, 0, 3,
                                                                       30000, 8634, 30006, 8877,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 30504, 0, 3,
                                                                       30006, 8637, 30012, 8886,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 30522, 0, 3,
                                                                       30012, 8640, 30018, 8895,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 30540, 0, 3,
                                                                       30018, 8643, 30024, 8904,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 30558, 0, 3,
                                                                       30024, 8646, 30030, 8913,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 30576, 0, 3,
                                                                       30030, 8649, 30036, 8922,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 30594, 0, 3,
                                                                       30036, 8652, 30042, 8931,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 30612, 0, 3,
                                                                       30042, 8655, 30048, 8940,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 30630, 0, 3,
                                                                       30048, 8658, 30054, 8949,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 30648, 0, 3,
                                                                       30054, 8661, 30060, 8958,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 30666, 0, 3,
                                                                       30060, 8664, 30066, 8967,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 30684, 0, 3,
                                                                       30072, 8670, 30090, 162,
                                                                       168, 8976, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 30720, 0, 3,
                                                                       30090, 8679, 30108, 168,
                                                                       174, 8994, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 30756, 0, 3,
                                                                       30108, 8688, 30126, 174,
                                                                       180, 9012, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 30792, 0, 3,
                                                                       30126, 8697, 30144, 180,
                                                                       186, 9030, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 30828, 0, 3,
                                                                       30144, 8706, 30162, 186,
                                                                       192, 9048, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 30864, 0, 3,
                                                                       30162, 8715, 30180, 192,
                                                                       198, 9066, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 30900, 0, 3,
                                                                       30180, 8724, 30198, 198,
                                                                       204, 9084, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 30936, 0, 3,
                                                                       30198, 8733, 30216, 204,
                                                                       210, 9102, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 30972, 0, 3,
                                                                       30216, 8742, 30234, 210,
                                                                       216, 9120, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 31008, 0, 3,
                                                                       30234, 8751, 30252, 216,
                                                                       222, 9138, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 31044, 0, 3,
                                                                       30252, 8760, 30270, 222,
                                                                       228, 9156, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 31080, 0, 3,
                                                                       30270, 8769, 30288, 228,
                                                                       234, 9174, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 31116, 0, 3,
                                                                       30288, 8778, 30306, 234,
                                                                       240, 9192, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 31152, 0, 3,
                                                                       30306, 8787, 30324, 240,
                                                                       246, 9210, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 31188, 0, 3,
                                                                       30324, 8796, 30342, 246,
                                                                       252, 9228, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 31224, 0, 3,
                                                                       30342, 8805, 30360, 252,
                                                                       258, 9246, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 31260, 0, 3,
                                                                       30378, 8823, 30396, 270,
                                                                       276, 9264, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 31296, 0, 3,
                                                                       30396, 8832, 30414, 276,
                                                                       282, 9282, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 31332, 0, 3,
                                                                       30414, 8841, 30432, 282,
                                                                       288, 9300, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 31368, 0, 3,
                                                                       30432, 8850, 30450, 288,
                                                                       294, 9318, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 31404, 0, 3,
                                                                       30450, 8859, 30468, 294,
                                                                       300, 9336, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 31440, 0, 3,
                                                                       30468, 8868, 30486, 300,
                                                                       306, 9354, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 31476, 0, 3,
                                                                       30486, 8877, 30504, 306,
                                                                       312, 9372, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 31512, 0, 3,
                                                                       30504, 8886, 30522, 312,
                                                                       318, 9390, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 31548, 0, 3,
                                                                       30522, 8895, 30540, 318,
                                                                       324, 9408, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 31584, 0, 3,
                                                                       30540, 8904, 30558, 324,
                                                                       330, 9426, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 31620, 0, 3,
                                                                       30558, 8913, 30576, 330,
                                                                       336, 9444, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 31656, 0, 3,
                                                                       30576, 8922, 30594, 336,
                                                                       342, 9462, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 31692, 0, 3,
                                                                       30594, 8931, 30612, 342,
                                                                       348, 9480, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 31728, 0, 3,
                                                                       30612, 8940, 30630, 348,
                                                                       354, 9498, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 31764, 0, 3,
                                                                       30630, 8949, 30648, 354,
                                                                       360, 9516, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 31800, 0, 3,
                                                                       30648, 8958, 30666, 360,
                                                                       366, 9534, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 31836, 0, 3,
                                                                       30684, 8976, 30720, 378,
                                                                       388, 9552, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 31896, 0, 3,
                                                                       30720, 8994, 30756, 388,
                                                                       398, 9582, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 31956, 0, 3,
                                                                       30756, 9012, 30792, 398,
                                                                       408, 9612, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 32016, 0, 3,
                                                                       30792, 9030, 30828, 408,
                                                                       418, 9642, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 32076, 0, 3,
                                                                       30828, 9048, 30864, 418,
                                                                       428, 9672, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 32136, 0, 3,
                                                                       30864, 9066, 30900, 428,
                                                                       438, 9702, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 32196, 0, 3,
                                                                       30900, 9084, 30936, 438,
                                                                       448, 9732, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 32256, 0, 3,
                                                                       30936, 9102, 30972, 448,
                                                                       458, 9762, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 32316, 0, 3,
                                                                       30972, 9120, 31008, 458,
                                                                       468, 9792, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 32376, 0, 3,
                                                                       31008, 9138, 31044, 468,
                                                                       478, 9822, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 32436, 0, 3,
                                                                       31044, 9156, 31080, 478,
                                                                       488, 9852, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 32496, 0, 3,
                                                                       31080, 9174, 31116, 488,
                                                                       498, 9882, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 32556, 0, 3,
                                                                       31116, 9192, 31152, 498,
                                                                       508, 9912, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 32616, 0, 3,
                                                                       31152, 9210, 31188, 508,
                                                                       518, 9942, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 32676, 0, 3,
                                                                       31188, 9228, 31224, 518,
                                                                       528, 9972, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 32736, 0, 3,
                                                                       31260, 9264, 31296, 548,
                                                                       558, 10002, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 32796, 0, 3,
                                                                       31296, 9282, 31332, 558,
                                                                       568, 10032, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 32856, 0, 3,
                                                                       31332, 9300, 31368, 568,
                                                                       578, 10062, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 32916, 0, 3,
                                                                       31368, 9318, 31404, 578,
                                                                       588, 10092, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 32976, 0, 3,
                                                                       31404, 9336, 31440, 588,
                                                                       598, 10122, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 33036, 0, 3,
                                                                       31440, 9354, 31476, 598,
                                                                       608, 10152, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 33096, 0, 3,
                                                                       31476, 9372, 31512, 608,
                                                                       618, 10182, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 33156, 0, 3,
                                                                       31512, 9390, 31548, 618,
                                                                       628, 10212, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 33216, 0, 3,
                                                                       31548, 9408, 31584, 628,
                                                                       638, 10242, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 33276, 0, 3,
                                                                       31584, 9426, 31620, 638,
                                                                       648, 10272, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 33336, 0, 3,
                                                                       31620, 9444, 31656, 648,
                                                                       658, 10302, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 33396, 0, 3,
                                                                       31656, 9462, 31692, 658,
                                                                       668, 10332, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 33456, 0, 3,
                                                                       31692, 9480, 31728, 668,
                                                                       678, 10362, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 33516, 0, 3,
                                                                       31728, 9498, 31764, 678,
                                                                       688, 10392, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 33576, 0, 3,
                                                                       31764, 9516, 31800, 688,
                                                                       698, 10422, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 33636, 0, 3,
                                                                       31836, 9552, 31896, 718,
                                                                       733, 10452, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 33726, 0, 3,
                                                                       31896, 9582, 31956, 733,
                                                                       748, 10497, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 33816, 0, 3,
                                                                       31956, 9612, 32016, 748,
                                                                       763, 10542, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 33906, 0, 3,
                                                                       32016, 9642, 32076, 763,
                                                                       778, 10587, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 33996, 0, 3,
                                                                       32076, 9672, 32136, 778,
                                                                       793, 10632, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 34086, 0, 3,
                                                                       32136, 9702, 32196, 793,
                                                                       808, 10677, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 34176, 0, 3,
                                                                       32196, 9732, 32256, 808,
                                                                       823, 10722, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 34266, 0, 3,
                                                                       32256, 9762, 32316, 823,
                                                                       838, 10767, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 34356, 0, 3,
                                                                       32316, 9792, 32376, 838,
                                                                       853, 10812, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 34446, 0, 3,
                                                                       32376, 9822, 32436, 853,
                                                                       868, 10857, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 34536, 0, 3,
                                                                       32436, 9852, 32496, 868,
                                                                       883, 10902, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 34626, 0, 3,
                                                                       32496, 9882, 32556, 883,
                                                                       898, 10947, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 34716, 0, 3,
                                                                       32556, 9912, 32616, 898,
                                                                       913, 10992, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 34806, 0, 3,
                                                                       32616, 9942, 32676, 913,
                                                                       928, 11037, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 34896, 0, 3,
                                                                       32736, 10002, 32796, 958,
                                                                       973, 11082, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 34986, 0, 3,
                                                                       32796, 10032, 32856, 973,
                                                                       988, 11127, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 35076, 0, 3,
                                                                       32856, 10062, 32916, 988,
                                                                       1003, 11172, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 35166, 0, 3,
                                                                       32916, 10092, 32976, 1003,
                                                                       1018, 11217, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 35256, 0, 3,
                                                                       32976, 10122, 33036, 1018,
                                                                       1033, 11262, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 35346, 0, 3,
                                                                       33036, 10152, 33096, 1033,
                                                                       1048, 11307, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 35436, 0, 3,
                                                                       33096, 10182, 33156, 1048,
                                                                       1063, 11352, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 35526, 0, 3,
                                                                       33156, 10212, 33216, 1063,
                                                                       1078, 11397, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 35616, 0, 3,
                                                                       33216, 10242, 33276, 1078,
                                                                       1093, 11442, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 35706, 0, 3,
                                                                       33276, 10272, 33336, 1093,
                                                                       1108, 11487, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 35796, 0, 3,
                                                                       33336, 10302, 33396, 1108,
                                                                       1123, 11532, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 35886, 0, 3,
                                                                       33396, 10332, 33456, 1123,
                                                                       1138, 11577, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 35976, 0, 3,
                                                                       33456, 10362, 33516, 1138,
                                                                       1153, 11622, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 36066, 0, 3,
                                                                       33516, 10392, 33576, 1153,
                                                                       1168, 11667, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 36156, 0, 3,
                                                                       33636, 10452, 33726, 1198,
                                                                       1219, 11712, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 36282, 0, 3,
                                                                       33726, 10497, 33816, 1219,
                                                                       1240, 11775, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 36408, 0, 3,
                                                                       33816, 10542, 33906, 1240,
                                                                       1261, 11838, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 36534, 0, 3,
                                                                       33906, 10587, 33996, 1261,
                                                                       1282, 11901, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 36660, 0, 3,
                                                                       33996, 10632, 34086, 1282,
                                                                       1303, 11964, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 36786, 0, 3,
                                                                       34086, 10677, 34176, 1303,
                                                                       1324, 12027, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 36912, 0, 3,
                                                                       34176, 10722, 34266, 1324,
                                                                       1345, 12090, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 37038, 0, 3,
                                                                       34266, 10767, 34356, 1345,
                                                                       1366, 12153, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 37164, 0, 3,
                                                                       34356, 10812, 34446, 1366,
                                                                       1387, 12216, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 37290, 0, 3,
                                                                       34446, 10857, 34536, 1387,
                                                                       1408, 12279, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 37416, 0, 3,
                                                                       34536, 10902, 34626, 1408,
                                                                       1429, 12342, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 37542, 0, 3,
                                                                       34626, 10947, 34716, 1429,
                                                                       1450, 12405, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 37668, 0, 3,
                                                                       34716, 10992, 34806, 1450,
                                                                       1471, 12468, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 37794, 0, 3,
                                                                       34896, 11082, 34986, 1513,
                                                                       1534, 12531, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 37920, 0, 3,
                                                                       34986, 11127, 35076, 1534,
                                                                       1555, 12594, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 38046, 0, 3,
                                                                       35076, 11172, 35166, 1555,
                                                                       1576, 12657, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 38172, 0, 3,
                                                                       35166, 11217, 35256, 1576,
                                                                       1597, 12720, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 38298, 0, 3,
                                                                       35256, 11262, 35346, 1597,
                                                                       1618, 12783, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 38424, 0, 3,
                                                                       35346, 11307, 35436, 1618,
                                                                       1639, 12846, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 38550, 0, 3,
                                                                       35436, 11352, 35526, 1639,
                                                                       1660, 12909, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 38676, 0, 3,
                                                                       35526, 11397, 35616, 1660,
                                                                       1681, 12972, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 38802, 0, 3,
                                                                       35616, 11442, 35706, 1681,
                                                                       1702, 13035, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 38928, 0, 3,
                                                                       35706, 11487, 35796, 1702,
                                                                       1723, 13098, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 39054, 0, 3,
                                                                       35796, 11532, 35886, 1723,
                                                                       1744, 13161, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 39180, 0, 3,
                                                                       35886, 11577, 35976, 1744,
                                                                       1765, 13224, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 39306, 0, 3,
                                                                       35976, 11622, 36066, 1765,
                                                                       1786, 13287, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 39432, 0, 3,
                                                                       36156, 11712, 36282, 1828,
                                                                       1856, 13350, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 39600, 0, 3,
                                                                       36282, 11775, 36408, 1856,
                                                                       1884, 13434, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 39768, 0, 3,
                                                                       36408, 11838, 36534, 1884,
                                                                       1912, 13518, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 39936, 0, 3,
                                                                       36534, 11901, 36660, 1912,
                                                                       1940, 13602, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 40104, 0, 3,
                                                                       36660, 11964, 36786, 1940,
                                                                       1968, 13686, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 40272, 0, 3,
                                                                       36786, 12027, 36912, 1968,
                                                                       1996, 13770, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 40440, 0, 3,
                                                                       36912, 12090, 37038, 1996,
                                                                       2024, 13854, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 40608, 0, 3,
                                                                       37038, 12153, 37164, 2024,
                                                                       2052, 13938, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 40776, 0, 3,
                                                                       37164, 12216, 37290, 2052,
                                                                       2080, 14022, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 40944, 0, 3,
                                                                       37290, 12279, 37416, 2080,
                                                                       2108, 14106, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 41112, 0, 3,
                                                                       37416, 12342, 37542, 2108,
                                                                       2136, 14190, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 41280, 0, 3,
                                                                       37542, 12405, 37668, 2136,
                                                                       2164, 14274, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 41448, 0, 3,
                                                                       37794, 12531, 37920, 2220,
                                                                       2248, 14358, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 41616, 0, 3,
                                                                       37920, 12594, 38046, 2248,
                                                                       2276, 14442, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 41784, 0, 3,
                                                                       38046, 12657, 38172, 2276,
                                                                       2304, 14526, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 41952, 0, 3,
                                                                       38172, 12720, 38298, 2304,
                                                                       2332, 14610, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 42120, 0, 3,
                                                                       38298, 12783, 38424, 2332,
                                                                       2360, 14694, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 42288, 0, 3,
                                                                       38424, 12846, 38550, 2360,
                                                                       2388, 14778, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 42456, 0, 3,
                                                                       38550, 12909, 38676, 2388,
                                                                       2416, 14862, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 42624, 0, 3,
                                                                       38676, 12972, 38802, 2416,
                                                                       2444, 14946, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 42792, 0, 3,
                                                                       38802, 13035, 38928, 2444,
                                                                       2472, 15030, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 42960, 0, 3,
                                                                       38928, 13098, 39054, 2472,
                                                                       2500, 15114, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 43128, 0, 3,
                                                                       39054, 13161, 39180, 2500,
                                                                       2528, 15198, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 43296, 0, 3,
                                                                       39180, 13224, 39306, 2528,
                                                                       2556, 15282, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 43464, 0, 3,
                                                                       39432, 13350, 39600, 2612,
                                                                       2648, 15366, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 43680, 0, 3,
                                                                       39600, 13434, 39768, 2648,
                                                                       2684, 15474, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 43896, 0, 3,
                                                                       39768, 13518, 39936, 2684,
                                                                       2720, 15582, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 44112, 0, 3,
                                                                       39936, 13602, 40104, 2720,
                                                                       2756, 15690, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 44328, 0, 3,
                                                                       40104, 13686, 40272, 2756,
                                                                       2792, 15798, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 44544, 0, 3,
                                                                       40272, 13770, 40440, 2792,
                                                                       2828, 15906, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 44760, 0, 3,
                                                                       40440, 13854, 40608, 2828,
                                                                       2864, 16014, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 44976, 0, 3,
                                                                       40608, 13938, 40776, 2864,
                                                                       2900, 16122, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 45192, 0, 3,
                                                                       40776, 14022, 40944, 2900,
                                                                       2936, 16230, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 45408, 0, 3,
                                                                       40944, 14106, 41112, 2936,
                                                                       2972, 16338, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 45624, 0, 3,
                                                                       41112, 14190, 41280, 2972,
                                                                       3008, 16446, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 45840, 0, 3,
                                                                       41448, 14358, 41616, 3080,
                                                                       3116, 16554, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 46056, 0, 3,
                                                                       41616, 14442, 41784, 3116,
                                                                       3152, 16662, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 46272, 0, 3,
                                                                       41784, 14526, 41952, 3152,
                                                                       3188, 16770, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 46488, 0, 3,
                                                                       41952, 14610, 42120, 3188,
                                                                       3224, 16878, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 46704, 0, 3,
                                                                       42120, 14694, 42288, 3224,
                                                                       3260, 16986, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 46920, 0, 3,
                                                                       42288, 14778, 42456, 3260,
                                                                       3296, 17094, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 47136, 0, 3,
                                                                       42456, 14862, 42624, 3296,
                                                                       3332, 17202, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 47352, 0, 3,
                                                                       42624, 14946, 42792, 3332,
                                                                       3368, 17310, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 47568, 0, 3,
                                                                       42792, 15030, 42960, 3368,
                                                                       3404, 17418, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 47784, 0, 3,
                                                                       42960, 15114, 43128, 3404,
                                                                       3440, 17526, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 48000, 0, 3,
                                                                       43128, 15198, 43296, 3440,
                                                                       3476, 17634, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 48216, 0, 3,
                                                                       43464, 15366, 43680, 3548,
                                                                       3593, 17742, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 48486, 0, 3,
                                                                       43680, 15474, 43896, 3593,
                                                                       3638, 17877, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 48756, 0, 3,
                                                                       43896, 15582, 44112, 3638,
                                                                       3683, 18012, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 49026, 0, 3,
                                                                       44112, 15690, 44328, 3683,
                                                                       3728, 18147, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 49296, 0, 3,
                                                                       44328, 15798, 44544, 3728,
                                                                       3773, 18282, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 49566, 0, 3,
                                                                       44544, 15906, 44760, 3773,
                                                                       3818, 18417, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 49836, 0, 3,
                                                                       44760, 16014, 44976, 3818,
                                                                       3863, 18552, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 50106, 0, 3,
                                                                       44976, 16122, 45192, 3863,
                                                                       3908, 18687, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 50376, 0, 3,
                                                                       45192, 16230, 45408, 3908,
                                                                       3953, 18822, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 50646, 0, 3,
                                                                       45408, 16338, 45624, 3953,
                                                                       3998, 18957, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 50916, 0, 3,
                                                                       45840, 16554, 46056, 4088,
                                                                       4133, 19092, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 51186, 0, 3,
                                                                       46056, 16662, 46272, 4133,
                                                                       4178, 19227, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 51456, 0, 3,
                                                                       46272, 16770, 46488, 4178,
                                                                       4223, 19362, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 51726, 0, 3,
                                                                       46488, 16878, 46704, 4223,
                                                                       4268, 19497, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 51996, 0, 3,
                                                                       46704, 16986, 46920, 4268,
                                                                       4313, 19632, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 52266, 0, 3,
                                                                       46920, 17094, 47136, 4313,
                                                                       4358, 19767, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 52536, 0, 3,
                                                                       47136, 17202, 47352, 4358,
                                                                       4403, 19902, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 52806, 0, 3,
                                                                       47352, 17310, 47568, 4403,
                                                                       4448, 20037, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 53076, 0, 3,
                                                                       47568, 17418, 47784, 4448,
                                                                       4493, 20172, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 53346, 0, 3,
                                                                       47784, 17526, 48000, 4493,
                                                                       4538, 20307, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 53616, 0, 3,
                                                                       48216, 17742, 48486, 4628,
                                                                       4683, 20442, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 53946, 0, 3,
                                                                       48486, 17877, 48756, 4683,
                                                                       4738, 20607, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 54276, 0, 3,
                                                                       48756, 18012, 49026, 4738,
                                                                       4793, 20772, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 54606, 0, 3,
                                                                       49026, 18147, 49296, 4793,
                                                                       4848, 20937, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 54936, 0, 3,
                                                                       49296, 18282, 49566, 4848,
                                                                       4903, 21102, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 55266, 0, 3,
                                                                       49566, 18417, 49836, 4903,
                                                                       4958, 21267, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 55596, 0, 3,
                                                                       49836, 18552, 50106, 4958,
                                                                       5013, 21432, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 55926, 0, 3,
                                                                       50106, 18687, 50376, 5013,
                                                                       5068, 21597, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 56256, 0, 3,
                                                                       50376, 18822, 50646, 5068,
                                                                       5123, 21762, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 56586, 0, 3,
                                                                       50916, 19092, 51186, 5233,
                                                                       5288, 21927, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 56916, 0, 3,
                                                                       51186, 19227, 51456, 5288,
                                                                       5343, 22092, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 57246, 0, 3,
                                                                       51456, 19362, 51726, 5343,
                                                                       5398, 22257, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 57576, 0, 3,
                                                                       51726, 19497, 51996, 5398,
                                                                       5453, 22422, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 57906, 0, 3,
                                                                       51996, 19632, 52266, 5453,
                                                                       5508, 22587, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 58236, 0, 3,
                                                                       52266, 19767, 52536, 5508,
                                                                       5563, 22752, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 58566, 0, 3,
                                                                       52536, 19902, 52806, 5563,
                                                                       5618, 22917, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 58896, 0, 3,
                                                                       52806, 20037, 53076, 5618,
                                                                       5673, 23082, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 59226, 0, 3,
                                                                       53076, 20172, 53346, 5673,
                                                                       5728, 23247, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 59556, 0, 3,
                                                                       53616, 20442, 53946, 5838,
                                                                       5904, 23412, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 59952, 0, 3,
                                                                       53946, 20607, 54276, 5904,
                                                                       5970, 23610, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 60348, 0, 3,
                                                                       54276, 20772, 54606, 5970,
                                                                       6036, 23808, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 60744, 0, 3,
                                                                       54606, 20937, 54936, 6036,
                                                                       6102, 24006, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 61140, 0, 3,
                                                                       54936, 21102, 55266, 6102,
                                                                       6168, 24204, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 61536, 0, 3,
                                                                       55266, 21267, 55596, 6168,
                                                                       6234, 24402, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 61932, 0, 3,
                                                                       55596, 21432, 55926, 6234,
                                                                       6300, 24600, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 62328, 0, 3,
                                                                       55926, 21597, 56256, 6300,
                                                                       6366, 24798, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 62724, 0, 3,
                                                                       56586, 21927, 56916, 6498,
                                                                       6564, 24996, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 63120, 0, 3,
                                                                       56916, 22092, 57246, 6564,
                                                                       6630, 25194, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 63516, 0, 3,
                                                                       57246, 22257, 57576, 6630,
                                                                       6696, 25392, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 63912, 0, 3,
                                                                       57576, 22422, 57906, 6696,
                                                                       6762, 25590, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 64308, 0, 3,
                                                                       57906, 22587, 58236, 6762,
                                                                       6828, 25788, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 64704, 0, 3,
                                                                       58236, 22752, 58566, 6828,
                                                                       6894, 25986, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 65100, 0, 3,
                                                                       58566, 22917, 58896, 6894,
                                                                       6960, 26184, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 65496, 0, 3,
                                                                       58896, 23082, 59226, 6960,
                                                                       7026, 26382, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 65892, 0, 3,
                                                                       59556, 23412, 59952, 7158,
                                                                       7236, 26580, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 66360, 0, 3,
                                                                       59952, 23610, 60348, 7236,
                                                                       7314, 26814, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 66828, 0, 3,
                                                                       60348, 23808, 60744, 7314,
                                                                       7392, 27048, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 67296, 0, 3,
                                                                       60744, 24006, 61140, 7392,
                                                                       7470, 27282, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 67764, 0, 3,
                                                                       61140, 24204, 61536, 7470,
                                                                       7548, 27516, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 68232, 0, 3,
                                                                       61536, 24402, 61932, 7548,
                                                                       7626, 27750, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 68700, 0, 3,
                                                                       61932, 24600, 62328, 7626,
                                                                       7704, 27984, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 69168, 0, 3,
                                                                       62724, 24996, 63120, 7860,
                                                                       7938, 28218, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 69636, 0, 3,
                                                                       63120, 25194, 63516, 7938,
                                                                       8016, 28452, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 70104, 0, 3,
                                                                       63516, 25392, 63912, 8016,
                                                                       8094, 28686, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 70572, 0, 3,
                                                                       63912, 25590, 64308, 8094,
                                                                       8172, 28920, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 71040, 0, 3,
                                                                       64308, 25788, 64704, 8172,
                                                                       8250, 29154, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 71508, 0, 3,
                                                                       64704, 25986, 65100, 8250,
                                                                       8328, 29388, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 71976, 0, 3,
                                                                       65100, 26184, 65496, 8328,
                                                                       8406, 29622, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 72444, 3, 8562,
                                                                       8565, 29868, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 72454, 3, 8565,
                                                                       8568, 29874, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 72464, 3, 8568,
                                                                       8571, 29880, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 72474, 3, 8571,
                                                                       8574, 29886, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 72484, 3, 8574,
                                                                       8577, 29892, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 72494, 3, 8577,
                                                                       8580, 29898, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 72504, 3, 8580,
                                                                       8583, 29904, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 72514, 3, 8583,
                                                                       8586, 29910, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 72524, 3, 8586,
                                                                       8589, 29916, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 72534, 3, 8589,
                                                                       8592, 29922, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 72544, 3, 8592,
                                                                       8595, 29928, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 72554, 3, 8595,
                                                                       8598, 29934, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 72564, 3, 8598,
                                                                       8601, 29940, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 72574, 3, 8601,
                                                                       8604, 29946, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 72584, 3, 8604,
                                                                       8607, 29952, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 72594, 3, 8607,
                                                                       8610, 29958, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 72604, 3, 8616,
                                                                       8619, 29976, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 72614, 3, 8619,
                                                                       8622, 29982, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 72624, 3, 8622,
                                                                       8625, 29988, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 72634, 3, 8625,
                                                                       8628, 29994, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 72644, 3, 8628,
                                                                       8631, 30000, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 72654, 3, 8631,
                                                                       8634, 30006, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 72664, 3, 8634,
                                                                       8637, 30012, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 72674, 3, 8637,
                                                                       8640, 30018, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 72684, 3, 8640,
                                                                       8643, 30024, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 72694, 3, 8643,
                                                                       8646, 30030, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 72704, 3, 8646,
                                                                       8649, 30036, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 72714, 3, 8649,
                                                                       8652, 30042, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 72724, 3, 8652,
                                                                       8655, 30048, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 72734, 3, 8655,
                                                                       8658, 30054, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 72744, 3, 8658,
                                                                       8661, 30060, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 72754, 3, 8661,
                                                                       8664, 30066, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 72764, 0, 3,
                                                                       72444, 29868, 72454,
                                                                       30108, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 72794, 0, 3,
                                                                       72454, 29874, 72464,
                                                                       30126, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 72824, 0, 3,
                                                                       72464, 29880, 72474,
                                                                       30144, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 72854, 0, 3,
                                                                       72474, 29886, 72484,
                                                                       30162, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 72884, 0, 3,
                                                                       72484, 29892, 72494,
                                                                       30180, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 72914, 0, 3,
                                                                       72494, 29898, 72504,
                                                                       30198, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 72944, 0, 3,
                                                                       72504, 29904, 72514,
                                                                       30216, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 72974, 0, 3,
                                                                       72514, 29910, 72524,
                                                                       30234, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 73004, 0, 3,
                                                                       72524, 29916, 72534,
                                                                       30252, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 73034, 0, 3,
                                                                       72534, 29922, 72544,
                                                                       30270, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 73064, 0, 3,
                                                                       72544, 29928, 72554,
                                                                       30288, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 73094, 0, 3,
                                                                       72554, 29934, 72564,
                                                                       30306, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 73124, 0, 3,
                                                                       72564, 29940, 72574,
                                                                       30324, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 73154, 0, 3,
                                                                       72574, 29946, 72584,
                                                                       30342, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 73184, 0, 3,
                                                                       72584, 29952, 72594,
                                                                       30360, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 73214, 0, 3,
                                                                       72604, 29976, 72614,
                                                                       30414, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 73244, 0, 3,
                                                                       72614, 29982, 72624,
                                                                       30432, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 73274, 0, 3,
                                                                       72624, 29988, 72634,
                                                                       30450, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 73304, 0, 3,
                                                                       72634, 29994, 72644,
                                                                       30468, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 73334, 0, 3,
                                                                       72644, 30000, 72654,
                                                                       30486, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 73364, 0, 3,
                                                                       72654, 30006, 72664,
                                                                       30504, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 73394, 0, 3,
                                                                       72664, 30012, 72674,
                                                                       30522, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 73424, 0, 3,
                                                                       72674, 30018, 72684,
                                                                       30540, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 73454, 0, 3,
                                                                       72684, 30024, 72694,
                                                                       30558, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 73484, 0, 3,
                                                                       72694, 30030, 72704,
                                                                       30576, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 73514, 0, 3,
                                                                       72704, 30036, 72714,
                                                                       30594, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 73544, 0, 3,
                                                                       72714, 30042, 72724,
                                                                       30612, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 73574, 0, 3,
                                                                       72724, 30048, 72734,
                                                                       30630, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 73604, 0, 3,
                                                                       72734, 30054, 72744,
                                                                       30648, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 73634, 0, 3,
                                                                       72744, 30060, 72754,
                                                                       30666, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 73664, 0, 3,
                                                                       72764, 30108, 72794, 8976,
                                                                       8994, 30756, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 73724, 0, 3,
                                                                       72794, 30126, 72824, 8994,
                                                                       9012, 30792, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 73784, 0, 3,
                                                                       72824, 30144, 72854, 9012,
                                                                       9030, 30828, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 73844, 0, 3,
                                                                       72854, 30162, 72884, 9030,
                                                                       9048, 30864, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 73904, 0, 3,
                                                                       72884, 30180, 72914, 9048,
                                                                       9066, 30900, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 73964, 0, 3,
                                                                       72914, 30198, 72944, 9066,
                                                                       9084, 30936, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 74024, 0, 3,
                                                                       72944, 30216, 72974, 9084,
                                                                       9102, 30972, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 74084, 0, 3,
                                                                       72974, 30234, 73004, 9102,
                                                                       9120, 31008, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 74144, 0, 3,
                                                                       73004, 30252, 73034, 9120,
                                                                       9138, 31044, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 74204, 0, 3,
                                                                       73034, 30270, 73064, 9138,
                                                                       9156, 31080, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 74264, 0, 3,
                                                                       73064, 30288, 73094, 9156,
                                                                       9174, 31116, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 74324, 0, 3,
                                                                       73094, 30306, 73124, 9174,
                                                                       9192, 31152, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 74384, 0, 3,
                                                                       73124, 30324, 73154, 9192,
                                                                       9210, 31188, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 74444, 0, 3,
                                                                       73154, 30342, 73184, 9210,
                                                                       9228, 31224, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 74504, 0, 3,
                                                                       73214, 30414, 73244, 9264,
                                                                       9282, 31332, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 74564, 0, 3,
                                                                       73244, 30432, 73274, 9282,
                                                                       9300, 31368, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 74624, 0, 3,
                                                                       73274, 30450, 73304, 9300,
                                                                       9318, 31404, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 74684, 0, 3,
                                                                       73304, 30468, 73334, 9318,
                                                                       9336, 31440, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 74744, 0, 3,
                                                                       73334, 30486, 73364, 9336,
                                                                       9354, 31476, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 74804, 0, 3,
                                                                       73364, 30504, 73394, 9354,
                                                                       9372, 31512, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 74864, 0, 3,
                                                                       73394, 30522, 73424, 9372,
                                                                       9390, 31548, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 74924, 0, 3,
                                                                       73424, 30540, 73454, 9390,
                                                                       9408, 31584, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 74984, 0, 3,
                                                                       73454, 30558, 73484, 9408,
                                                                       9426, 31620, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 75044, 0, 3,
                                                                       73484, 30576, 73514, 9426,
                                                                       9444, 31656, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 75104, 0, 3,
                                                                       73514, 30594, 73544, 9444,
                                                                       9462, 31692, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 75164, 0, 3,
                                                                       73544, 30612, 73574, 9462,
                                                                       9480, 31728, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 75224, 0, 3,
                                                                       73574, 30630, 73604, 9480,
                                                                       9498, 31764, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 75284, 0, 3,
                                                                       73604, 30648, 73634, 9498,
                                                                       9516, 31800, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 75344, 0, 3,
                                                                       73664, 30756, 73724, 9552,
                                                                       9582, 31956, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 75444, 0, 3,
                                                                       73724, 30792, 73784, 9582,
                                                                       9612, 32016, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 75544, 0, 3,
                                                                       73784, 30828, 73844, 9612,
                                                                       9642, 32076, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 75644, 0, 3,
                                                                       73844, 30864, 73904, 9642,
                                                                       9672, 32136, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 75744, 0, 3,
                                                                       73904, 30900, 73964, 9672,
                                                                       9702, 32196, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 75844, 0, 3,
                                                                       73964, 30936, 74024, 9702,
                                                                       9732, 32256, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 75944, 0, 3,
                                                                       74024, 30972, 74084, 9732,
                                                                       9762, 32316, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 76044, 0, 3,
                                                                       74084, 31008, 74144, 9762,
                                                                       9792, 32376, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 76144, 0, 3,
                                                                       74144, 31044, 74204, 9792,
                                                                       9822, 32436, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 76244, 0, 3,
                                                                       74204, 31080, 74264, 9822,
                                                                       9852, 32496, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 76344, 0, 3,
                                                                       74264, 31116, 74324, 9852,
                                                                       9882, 32556, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 76444, 0, 3,
                                                                       74324, 31152, 74384, 9882,
                                                                       9912, 32616, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 76544, 0, 3,
                                                                       74384, 31188, 74444, 9912,
                                                                       9942, 32676, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 76644, 0, 3,
                                                                       74504, 31332, 74564,
                                                                       10002, 10032, 32856,
                                                                       ncols, gamma, p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 76744, 0, 3,
                                                                       74564, 31368, 74624,
                                                                       10032, 10062, 32916,
                                                                       ncols, gamma, p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 76844, 0, 3,
                                                                       74624, 31404, 74684,
                                                                       10062, 10092, 32976,
                                                                       ncols, gamma, p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 76944, 0, 3,
                                                                       74684, 31440, 74744,
                                                                       10092, 10122, 33036,
                                                                       ncols, gamma, p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 77044, 0, 3,
                                                                       74744, 31476, 74804,
                                                                       10122, 10152, 33096,
                                                                       ncols, gamma, p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 77144, 0, 3,
                                                                       74804, 31512, 74864,
                                                                       10152, 10182, 33156,
                                                                       ncols, gamma, p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 77244, 0, 3,
                                                                       74864, 31548, 74924,
                                                                       10182, 10212, 33216,
                                                                       ncols, gamma, p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 77344, 0, 3,
                                                                       74924, 31584, 74984,
                                                                       10212, 10242, 33276,
                                                                       ncols, gamma, p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 77444, 0, 3,
                                                                       74984, 31620, 75044,
                                                                       10242, 10272, 33336,
                                                                       ncols, gamma, p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 77544, 0, 3,
                                                                       75044, 31656, 75104,
                                                                       10272, 10302, 33396,
                                                                       ncols, gamma, p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 77644, 0, 3,
                                                                       75104, 31692, 75164,
                                                                       10302, 10332, 33456,
                                                                       ncols, gamma, p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 77744, 0, 3,
                                                                       75164, 31728, 75224,
                                                                       10332, 10362, 33516,
                                                                       ncols, gamma, p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 77844, 0, 3,
                                                                       75224, 31764, 75284,
                                                                       10362, 10392, 33576,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 77944, 0, 3,
                                                                       75344, 31956, 75444,
                                                                       10452, 10497, 33816,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 78094, 0, 3,
                                                                       75444, 32016, 75544,
                                                                       10497, 10542, 33906,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 78244, 0, 3,
                                                                       75544, 32076, 75644,
                                                                       10542, 10587, 33996,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 78394, 0, 3,
                                                                       75644, 32136, 75744,
                                                                       10587, 10632, 34086,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 78544, 0, 3,
                                                                       75744, 32196, 75844,
                                                                       10632, 10677, 34176,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 78694, 0, 3,
                                                                       75844, 32256, 75944,
                                                                       10677, 10722, 34266,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 78844, 0, 3,
                                                                       75944, 32316, 76044,
                                                                       10722, 10767, 34356,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 78994, 0, 3,
                                                                       76044, 32376, 76144,
                                                                       10767, 10812, 34446,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 79144, 0, 3,
                                                                       76144, 32436, 76244,
                                                                       10812, 10857, 34536,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 79294, 0, 3,
                                                                       76244, 32496, 76344,
                                                                       10857, 10902, 34626,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 79444, 0, 3,
                                                                       76344, 32556, 76444,
                                                                       10902, 10947, 34716,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 79594, 0, 3,
                                                                       76444, 32616, 76544,
                                                                       10947, 10992, 34806,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 79744, 0, 3,
                                                                       76644, 32856, 76744,
                                                                       11082, 11127, 35076,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 79894, 0, 3,
                                                                       76744, 32916, 76844,
                                                                       11127, 11172, 35166,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 80044, 0, 3,
                                                                       76844, 32976, 76944,
                                                                       11172, 11217, 35256,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 80194, 0, 3,
                                                                       76944, 33036, 77044,
                                                                       11217, 11262, 35346,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 80344, 0, 3,
                                                                       77044, 33096, 77144,
                                                                       11262, 11307, 35436,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 80494, 0, 3,
                                                                       77144, 33156, 77244,
                                                                       11307, 11352, 35526,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 80644, 0, 3,
                                                                       77244, 33216, 77344,
                                                                       11352, 11397, 35616,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 80794, 0, 3,
                                                                       77344, 33276, 77444,
                                                                       11397, 11442, 35706,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 80944, 0, 3,
                                                                       77444, 33336, 77544,
                                                                       11442, 11487, 35796,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 81094, 0, 3,
                                                                       77544, 33396, 77644,
                                                                       11487, 11532, 35886,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 81244, 0, 3,
                                                                       77644, 33456, 77744,
                                                                       11532, 11577, 35976,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 81394, 0, 3,
                                                                       77744, 33516, 77844,
                                                                       11577, 11622, 36066,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 81544, 0, 3,
                                                                       77944, 33816, 78094,
                                                                       11712, 11775, 36408,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 81754, 0, 3,
                                                                       78094, 33906, 78244,
                                                                       11775, 11838, 36534,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 81964, 0, 3,
                                                                       78244, 33996, 78394,
                                                                       11838, 11901, 36660,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 82174, 0, 3,
                                                                       78394, 34086, 78544,
                                                                       11901, 11964, 36786,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 82384, 0, 3,
                                                                       78544, 34176, 78694,
                                                                       11964, 12027, 36912,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 82594, 0, 3,
                                                                       78694, 34266, 78844,
                                                                       12027, 12090, 37038,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 82804, 0, 3,
                                                                       78844, 34356, 78994,
                                                                       12090, 12153, 37164,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 83014, 0, 3,
                                                                       78994, 34446, 79144,
                                                                       12153, 12216, 37290,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 83224, 0, 3,
                                                                       79144, 34536, 79294,
                                                                       12216, 12279, 37416,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 83434, 0, 3,
                                                                       79294, 34626, 79444,
                                                                       12279, 12342, 37542,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 83644, 0, 3,
                                                                       79444, 34716, 79594,
                                                                       12342, 12405, 37668,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 83854, 0, 3,
                                                                       79744, 35076, 79894,
                                                                       12531, 12594, 38046,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 84064, 0, 3,
                                                                       79894, 35166, 80044,
                                                                       12594, 12657, 38172,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 84274, 0, 3,
                                                                       80044, 35256, 80194,
                                                                       12657, 12720, 38298,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 84484, 0, 3,
                                                                       80194, 35346, 80344,
                                                                       12720, 12783, 38424,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 84694, 0, 3,
                                                                       80344, 35436, 80494,
                                                                       12783, 12846, 38550,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 84904, 0, 3,
                                                                       80494, 35526, 80644,
                                                                       12846, 12909, 38676,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 85114, 0, 3,
                                                                       80644, 35616, 80794,
                                                                       12909, 12972, 38802,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 85324, 0, 3,
                                                                       80794, 35706, 80944,
                                                                       12972, 13035, 38928,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 85534, 0, 3,
                                                                       80944, 35796, 81094,
                                                                       13035, 13098, 39054,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 85744, 0, 3,
                                                                       81094, 35886, 81244,
                                                                       13098, 13161, 39180,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 85954, 0, 3,
                                                                       81244, 35976, 81394,
                                                                       13161, 13224, 39306,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 86164, 0, 3,
                                                                       81544, 36408, 81754,
                                                                       13350, 13434, 39768,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 86444, 0, 3,
                                                                       81754, 36534, 81964,
                                                                       13434, 13518, 39936,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 86724, 0, 3,
                                                                       81964, 36660, 82174,
                                                                       13518, 13602, 40104,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 87004, 0, 3,
                                                                       82174, 36786, 82384,
                                                                       13602, 13686, 40272,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 87284, 0, 3,
                                                                       82384, 36912, 82594,
                                                                       13686, 13770, 40440,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 87564, 0, 3,
                                                                       82594, 37038, 82804,
                                                                       13770, 13854, 40608,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 87844, 0, 3,
                                                                       82804, 37164, 83014,
                                                                       13854, 13938, 40776,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 88124, 0, 3,
                                                                       83014, 37290, 83224,
                                                                       13938, 14022, 40944,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 88404, 0, 3,
                                                                       83224, 37416, 83434,
                                                                       14022, 14106, 41112,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 88684, 0, 3,
                                                                       83434, 37542, 83644,
                                                                       14106, 14190, 41280,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 88964, 0, 3,
                                                                       83854, 38046, 84064,
                                                                       14358, 14442, 41784,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 89244, 0, 3,
                                                                       84064, 38172, 84274,
                                                                       14442, 14526, 41952,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 89524, 0, 3,
                                                                       84274, 38298, 84484,
                                                                       14526, 14610, 42120,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 89804, 0, 3,
                                                                       84484, 38424, 84694,
                                                                       14610, 14694, 42288,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 90084, 0, 3,
                                                                       84694, 38550, 84904,
                                                                       14694, 14778, 42456,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 90364, 0, 3,
                                                                       84904, 38676, 85114,
                                                                       14778, 14862, 42624,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 90644, 0, 3,
                                                                       85114, 38802, 85324,
                                                                       14862, 14946, 42792,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 90924, 0, 3,
                                                                       85324, 38928, 85534,
                                                                       14946, 15030, 42960,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 91204, 0, 3,
                                                                       85534, 39054, 85744,
                                                                       15030, 15114, 43128,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 91484, 0, 3,
                                                                       85744, 39180, 85954,
                                                                       15114, 15198, 43296,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 91764, 0, 3,
                                                                       86164, 39768, 86444,
                                                                       15366, 15474, 43896,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 92124, 0, 3,
                                                                       86444, 39936, 86724,
                                                                       15474, 15582, 44112,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 92484, 0, 3,
                                                                       86724, 40104, 87004,
                                                                       15582, 15690, 44328,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 92844, 0, 3,
                                                                       87004, 40272, 87284,
                                                                       15690, 15798, 44544,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 93204, 0, 3,
                                                                       87284, 40440, 87564,
                                                                       15798, 15906, 44760,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 93564, 0, 3,
                                                                       87564, 40608, 87844,
                                                                       15906, 16014, 44976,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 93924, 0, 3,
                                                                       87844, 40776, 88124,
                                                                       16014, 16122, 45192,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 94284, 0, 3,
                                                                       88124, 40944, 88404,
                                                                       16122, 16230, 45408,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 94644, 0, 3,
                                                                       88404, 41112, 88684,
                                                                       16230, 16338, 45624,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 95004, 0, 3,
                                                                       88964, 41784, 89244,
                                                                       16554, 16662, 46272,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 95364, 0, 3,
                                                                       89244, 41952, 89524,
                                                                       16662, 16770, 46488,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 95724, 0, 3,
                                                                       89524, 42120, 89804,
                                                                       16770, 16878, 46704,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 96084, 0, 3,
                                                                       89804, 42288, 90084,
                                                                       16878, 16986, 46920,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 96444, 0, 3,
                                                                       90084, 42456, 90364,
                                                                       16986, 17094, 47136,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 96804, 0, 3,
                                                                       90364, 42624, 90644,
                                                                       17094, 17202, 47352,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 97164, 0, 3,
                                                                       90644, 42792, 90924,
                                                                       17202, 17310, 47568,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 97524, 0, 3,
                                                                       90924, 42960, 91204,
                                                                       17310, 17418, 47784,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 97884, 0, 3,
                                                                       91204, 43128, 91484,
                                                                       17418, 17526, 48000,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 98244, 0, 3,
                                                                       91764, 43896, 92124,
                                                                       17742, 17877, 48756,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 98694, 0, 3,
                                                                       92124, 44112, 92484,
                                                                       17877, 18012, 49026,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 99144, 0, 3,
                                                                       92484, 44328, 92844,
                                                                       18012, 18147, 49296,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 99594, 0, 3,
                                                                       92844, 44544, 93204,
                                                                       18147, 18282, 49566,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 100044, 0, 3,
                                                                       93204, 44760, 93564,
                                                                       18282, 18417, 49836,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 100494, 0, 3,
                                                                       93564, 44976, 93924,
                                                                       18417, 18552, 50106,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 100944, 0, 3,
                                                                       93924, 45192, 94284,
                                                                       18552, 18687, 50376,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 101394, 0, 3,
                                                                       94284, 45408, 94644,
                                                                       18687, 18822, 50646,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 101844, 0, 3,
                                                                       95004, 46272, 95364,
                                                                       19092, 19227, 51456,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 102294, 0, 3,
                                                                       95364, 46488, 95724,
                                                                       19227, 19362, 51726,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 102744, 0, 3,
                                                                       95724, 46704, 96084,
                                                                       19362, 19497, 51996,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 103194, 0, 3,
                                                                       96084, 46920, 96444,
                                                                       19497, 19632, 52266,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 103644, 0, 3,
                                                                       96444, 47136, 96804,
                                                                       19632, 19767, 52536,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 104094, 0, 3,
                                                                       96804, 47352, 97164,
                                                                       19767, 19902, 52806,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 104544, 0, 3,
                                                                       97164, 47568, 97524,
                                                                       19902, 20037, 53076,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 104994, 0, 3,
                                                                       97524, 47784, 97884,
                                                                       20037, 20172, 53346,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 105444, 0, 3,
                                                                       98244, 48756, 98694,
                                                                       20442, 20607, 54276,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 105994, 0, 3,
                                                                       98694, 49026, 99144,
                                                                       20607, 20772, 54606,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 106544, 0, 3,
                                                                       99144, 49296, 99594,
                                                                       20772, 20937, 54936,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 107094, 0, 3,
                                                                       99594, 49566, 100044,
                                                                       20937, 21102, 55266,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 107644, 0, 3,
                                                                       100044, 49836, 100494,
                                                                       21102, 21267, 55596,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 108194, 0, 3,
                                                                       100494, 50106, 100944,
                                                                       21267, 21432, 55926,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 108744, 0, 3,
                                                                       100944, 50376, 101394,
                                                                       21432, 21597, 56256,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 109294, 0, 3,
                                                                       101844, 51456, 102294,
                                                                       21927, 22092, 57246,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 109844, 0, 3,
                                                                       102294, 51726, 102744,
                                                                       22092, 22257, 57576,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 110394, 0, 3,
                                                                       102744, 51996, 103194,
                                                                       22257, 22422, 57906,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 110944, 0, 3,
                                                                       103194, 52266, 103644,
                                                                       22422, 22587, 58236,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 111494, 0, 3,
                                                                       103644, 52536, 104094,
                                                                       22587, 22752, 58566,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 112044, 0, 3,
                                                                       104094, 52806, 104544,
                                                                       22752, 22917, 58896,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 112594, 0, 3,
                                                                       104544, 53076, 104994,
                                                                       22917, 23082, 59226,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 113144, 0, 3,
                                                                       105444, 54276, 105994,
                                                                       23412, 23610, 60348,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 113804, 0, 3,
                                                                       105994, 54606, 106544,
                                                                       23610, 23808, 60744,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 114464, 0, 3,
                                                                       106544, 54936, 107094,
                                                                       23808, 24006, 61140,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 115124, 0, 3,
                                                                       107094, 55266, 107644,
                                                                       24006, 24204, 61536,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 115784, 0, 3,
                                                                       107644, 55596, 108194,
                                                                       24204, 24402, 61932,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 116444, 0, 3,
                                                                       108194, 55926, 108744,
                                                                       24402, 24600, 62328,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 117104, 0, 3,
                                                                       109294, 57246, 109844,
                                                                       24996, 25194, 63516,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 117764, 0, 3,
                                                                       109844, 57576, 110394,
                                                                       25194, 25392, 63912,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 118424, 0, 3,
                                                                       110394, 57906, 110944,
                                                                       25392, 25590, 64308,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 119084, 0, 3,
                                                                       110944, 58236, 111494,
                                                                       25590, 25788, 64704,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 119744, 0, 3,
                                                                       111494, 58566, 112044,
                                                                       25788, 25986, 65100,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 120404, 0, 3,
                                                                       112044, 58896, 112594,
                                                                       25986, 26184, 65496,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 121064, 0, 3,
                                                                       113144, 60348, 113804,
                                                                       26580, 26814, 66828,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 121844, 0, 3,
                                                                       113804, 60744, 114464,
                                                                       26814, 27048, 67296,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 122624, 0, 3,
                                                                       114464, 61140, 115124,
                                                                       27048, 27282, 67764,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 123404, 0, 3,
                                                                       115124, 61536, 115784,
                                                                       27282, 27516, 68232,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 124184, 0, 3,
                                                                       115784, 61932, 116444,
                                                                       27516, 27750, 68700,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 124964, 0, 3,
                                                                       117104, 63516, 117764,
                                                                       28218, 28452, 70104,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 125744, 0, 3,
                                                                       117764, 63912, 118424,
                                                                       28452, 28686, 70572,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 126524, 0, 3,
                                                                       118424, 64308, 119084,
                                                                       28686, 28920, 71040,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 127304, 0, 3,
                                                                       119084, 64704, 119744,
                                                                       28920, 29154, 71508,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 128084, 0, 3,
                                                                       119744, 65100, 120404,
                                                                       29154, 29388, 71976,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 128864, 3, 29856,
                                                                       29862, 72444, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 128879, 3, 29862,
                                                                       29868, 72454, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 128894, 3, 29868,
                                                                       29874, 72464, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 128909, 3, 29874,
                                                                       29880, 72474, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 128924, 3, 29880,
                                                                       29886, 72484, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 128939, 3, 29886,
                                                                       29892, 72494, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 128954, 3, 29892,
                                                                       29898, 72504, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 128969, 3, 29898,
                                                                       29904, 72514, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 128984, 3, 29904,
                                                                       29910, 72524, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 128999, 3, 29910,
                                                                       29916, 72534, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129014, 3, 29916,
                                                                       29922, 72544, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129029, 3, 29922,
                                                                       29928, 72554, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129044, 3, 29928,
                                                                       29934, 72564, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129059, 3, 29934,
                                                                       29940, 72574, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129074, 3, 29940,
                                                                       29946, 72584, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129089, 3, 29946,
                                                                       29952, 72594, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129104, 3, 29964,
                                                                       29970, 72604, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129119, 3, 29970,
                                                                       29976, 72614, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129134, 3, 29976,
                                                                       29982, 72624, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129149, 3, 29982,
                                                                       29988, 72634, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129164, 3, 29988,
                                                                       29994, 72644, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129179, 3, 29994,
                                                                       30000, 72654, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129194, 3, 30000,
                                                                       30006, 72664, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129209, 3, 30006,
                                                                       30012, 72674, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129224, 3, 30012,
                                                                       30018, 72684, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129239, 3, 30018,
                                                                       30024, 72694, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129254, 3, 30024,
                                                                       30030, 72704, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129269, 3, 30030,
                                                                       30036, 72714, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129284, 3, 30036,
                                                                       30042, 72724, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129299, 3, 30042,
                                                                       30048, 72734, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129314, 3, 30048,
                                                                       30054, 72744, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129329, 3, 30054,
                                                                       30060, 72754, ncols,
                                                                       gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 129344, 0, 3,
                                                                       128864, 72444, 128879,
                                                                       30072, 30090, 72764,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 129389, 0, 3,
                                                                       128879, 72454, 128894,
                                                                       30090, 30108, 72794,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 129434, 0, 3,
                                                                       128894, 72464, 128909,
                                                                       30108, 30126, 72824,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 129479, 0, 3,
                                                                       128909, 72474, 128924,
                                                                       30126, 30144, 72854,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 129524, 0, 3,
                                                                       128924, 72484, 128939,
                                                                       30144, 30162, 72884,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 129569, 0, 3,
                                                                       128939, 72494, 128954,
                                                                       30162, 30180, 72914,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 129614, 0, 3,
                                                                       128954, 72504, 128969,
                                                                       30180, 30198, 72944,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 129659, 0, 3,
                                                                       128969, 72514, 128984,
                                                                       30198, 30216, 72974,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 129704, 0, 3,
                                                                       128984, 72524, 128999,
                                                                       30216, 30234, 73004,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 129749, 0, 3,
                                                                       128999, 72534, 129014,
                                                                       30234, 30252, 73034,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 129794, 0, 3,
                                                                       129014, 72544, 129029,
                                                                       30252, 30270, 73064,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 129839, 0, 3,
                                                                       129029, 72554, 129044,
                                                                       30270, 30288, 73094,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 129884, 0, 3,
                                                                       129044, 72564, 129059,
                                                                       30288, 30306, 73124,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 129929, 0, 3,
                                                                       129059, 72574, 129074,
                                                                       30306, 30324, 73154,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 129974, 0, 3,
                                                                       129074, 72584, 129089,
                                                                       30324, 30342, 73184,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 130019, 0, 3,
                                                                       129104, 72604, 129119,
                                                                       30378, 30396, 73214,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 130064, 0, 3,
                                                                       129119, 72614, 129134,
                                                                       30396, 30414, 73244,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 130109, 0, 3,
                                                                       129134, 72624, 129149,
                                                                       30414, 30432, 73274,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 130154, 0, 3,
                                                                       129149, 72634, 129164,
                                                                       30432, 30450, 73304,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 130199, 0, 3,
                                                                       129164, 72644, 129179,
                                                                       30450, 30468, 73334,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 130244, 0, 3,
                                                                       129179, 72654, 129194,
                                                                       30468, 30486, 73364,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 130289, 0, 3,
                                                                       129194, 72664, 129209,
                                                                       30486, 30504, 73394,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 130334, 0, 3,
                                                                       129209, 72674, 129224,
                                                                       30504, 30522, 73424,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 130379, 0, 3,
                                                                       129224, 72684, 129239,
                                                                       30522, 30540, 73454,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 130424, 0, 3,
                                                                       129239, 72694, 129254,
                                                                       30540, 30558, 73484,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 130469, 0, 3,
                                                                       129254, 72704, 129269,
                                                                       30558, 30576, 73514,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 130514, 0, 3,
                                                                       129269, 72714, 129284,
                                                                       30576, 30594, 73544,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 130559, 0, 3,
                                                                       129284, 72724, 129299,
                                                                       30594, 30612, 73574,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 130604, 0, 3,
                                                                       129299, 72734, 129314,
                                                                       30612, 30630, 73604,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 130649, 0, 3,
                                                                       129314, 72744, 129329,
                                                                       30630, 30648, 73634,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 130694, 0, 3,
                                                                       129344, 72764, 129389,
                                                                       30684, 30720, 73664,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 130784, 0, 3,
                                                                       129389, 72794, 129434,
                                                                       30720, 30756, 73724,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 130874, 0, 3,
                                                                       129434, 72824, 129479,
                                                                       30756, 30792, 73784,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 130964, 0, 3,
                                                                       129479, 72854, 129524,
                                                                       30792, 30828, 73844,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 131054, 0, 3,
                                                                       129524, 72884, 129569,
                                                                       30828, 30864, 73904,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 131144, 0, 3,
                                                                       129569, 72914, 129614,
                                                                       30864, 30900, 73964,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 131234, 0, 3,
                                                                       129614, 72944, 129659,
                                                                       30900, 30936, 74024,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 131324, 0, 3,
                                                                       129659, 72974, 129704,
                                                                       30936, 30972, 74084,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 131414, 0, 3,
                                                                       129704, 73004, 129749,
                                                                       30972, 31008, 74144,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 131504, 0, 3,
                                                                       129749, 73034, 129794,
                                                                       31008, 31044, 74204,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 131594, 0, 3,
                                                                       129794, 73064, 129839,
                                                                       31044, 31080, 74264,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 131684, 0, 3,
                                                                       129839, 73094, 129884,
                                                                       31080, 31116, 74324,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 131774, 0, 3,
                                                                       129884, 73124, 129929,
                                                                       31116, 31152, 74384,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 131864, 0, 3,
                                                                       129929, 73154, 129974,
                                                                       31152, 31188, 74444,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 131954, 0, 3,
                                                                       130019, 73214, 130064,
                                                                       31260, 31296, 74504,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 132044, 0, 3,
                                                                       130064, 73244, 130109,
                                                                       31296, 31332, 74564,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 132134, 0, 3,
                                                                       130109, 73274, 130154,
                                                                       31332, 31368, 74624,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 132224, 0, 3,
                                                                       130154, 73304, 130199,
                                                                       31368, 31404, 74684,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 132314, 0, 3,
                                                                       130199, 73334, 130244,
                                                                       31404, 31440, 74744,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 132404, 0, 3,
                                                                       130244, 73364, 130289,
                                                                       31440, 31476, 74804,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 132494, 0, 3,
                                                                       130289, 73394, 130334,
                                                                       31476, 31512, 74864,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 132584, 0, 3,
                                                                       130334, 73424, 130379,
                                                                       31512, 31548, 74924,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 132674, 0, 3,
                                                                       130379, 73454, 130424,
                                                                       31548, 31584, 74984,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 132764, 0, 3,
                                                                       130424, 73484, 130469,
                                                                       31584, 31620, 75044,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 132854, 0, 3,
                                                                       130469, 73514, 130514,
                                                                       31620, 31656, 75104,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 132944, 0, 3,
                                                                       130514, 73544, 130559,
                                                                       31656, 31692, 75164,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 133034, 0, 3,
                                                                       130559, 73574, 130604,
                                                                       31692, 31728, 75224,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 133124, 0, 3,
                                                                       130604, 73604, 130649,
                                                                       31728, 31764, 75284,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 133214, 0, 3,
                                                                       130694, 73664, 130784,
                                                                       31836, 31896, 75344,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 133364, 0, 3,
                                                                       130784, 73724, 130874,
                                                                       31896, 31956, 75444,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 133514, 0, 3,
                                                                       130874, 73784, 130964,
                                                                       31956, 32016, 75544,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 133664, 0, 3,
                                                                       130964, 73844, 131054,
                                                                       32016, 32076, 75644,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 133814, 0, 3,
                                                                       131054, 73904, 131144,
                                                                       32076, 32136, 75744,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 133964, 0, 3,
                                                                       131144, 73964, 131234,
                                                                       32136, 32196, 75844,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 134114, 0, 3,
                                                                       131234, 74024, 131324,
                                                                       32196, 32256, 75944,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 134264, 0, 3,
                                                                       131324, 74084, 131414,
                                                                       32256, 32316, 76044,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 134414, 0, 3,
                                                                       131414, 74144, 131504,
                                                                       32316, 32376, 76144,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 134564, 0, 3,
                                                                       131504, 74204, 131594,
                                                                       32376, 32436, 76244,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 134714, 0, 3,
                                                                       131594, 74264, 131684,
                                                                       32436, 32496, 76344,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 134864, 0, 3,
                                                                       131684, 74324, 131774,
                                                                       32496, 32556, 76444,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 135014, 0, 3,
                                                                       131774, 74384, 131864,
                                                                       32556, 32616, 76544,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 135164, 0, 3,
                                                                       131954, 74504, 132044,
                                                                       32736, 32796, 76644,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 135314, 0, 3,
                                                                       132044, 74564, 132134,
                                                                       32796, 32856, 76744,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 135464, 0, 3,
                                                                       132134, 74624, 132224,
                                                                       32856, 32916, 76844,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 135614, 0, 3,
                                                                       132224, 74684, 132314,
                                                                       32916, 32976, 76944,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 135764, 0, 3,
                                                                       132314, 74744, 132404,
                                                                       32976, 33036, 77044,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 135914, 0, 3,
                                                                       132404, 74804, 132494,
                                                                       33036, 33096, 77144,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 136064, 0, 3,
                                                                       132494, 74864, 132584,
                                                                       33096, 33156, 77244,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 136214, 0, 3,
                                                                       132584, 74924, 132674,
                                                                       33156, 33216, 77344,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 136364, 0, 3,
                                                                       132674, 74984, 132764,
                                                                       33216, 33276, 77444,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 136514, 0, 3,
                                                                       132764, 75044, 132854,
                                                                       33276, 33336, 77544,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 136664, 0, 3,
                                                                       132854, 75104, 132944,
                                                                       33336, 33396, 77644,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 136814, 0, 3,
                                                                       132944, 75164, 133034,
                                                                       33396, 33456, 77744,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 136964, 0, 3,
                                                                       133034, 75224, 133124,
                                                                       33456, 33516, 77844,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 137114, 0, 3,
                                                                       133214, 75344, 133364,
                                                                       33636, 33726, 77944,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 137339, 0, 3,
                                                                       133364, 75444, 133514,
                                                                       33726, 33816, 78094,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 137564, 0, 3,
                                                                       133514, 75544, 133664,
                                                                       33816, 33906, 78244,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 137789, 0, 3,
                                                                       133664, 75644, 133814,
                                                                       33906, 33996, 78394,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 138014, 0, 3,
                                                                       133814, 75744, 133964,
                                                                       33996, 34086, 78544,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 138239, 0, 3,
                                                                       133964, 75844, 134114,
                                                                       34086, 34176, 78694,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 138464, 0, 3,
                                                                       134114, 75944, 134264,
                                                                       34176, 34266, 78844,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 138689, 0, 3,
                                                                       134264, 76044, 134414,
                                                                       34266, 34356, 78994,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 138914, 0, 3,
                                                                       134414, 76144, 134564,
                                                                       34356, 34446, 79144,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 139139, 0, 3,
                                                                       134564, 76244, 134714,
                                                                       34446, 34536, 79294,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 139364, 0, 3,
                                                                       134714, 76344, 134864,
                                                                       34536, 34626, 79444,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 139589, 0, 3,
                                                                       134864, 76444, 135014,
                                                                       34626, 34716, 79594,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 139814, 0, 3,
                                                                       135164, 76644, 135314,
                                                                       34896, 34986, 79744,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 140039, 0, 3,
                                                                       135314, 76744, 135464,
                                                                       34986, 35076, 79894,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 140264, 0, 3,
                                                                       135464, 76844, 135614,
                                                                       35076, 35166, 80044,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 140489, 0, 3,
                                                                       135614, 76944, 135764,
                                                                       35166, 35256, 80194,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 140714, 0, 3,
                                                                       135764, 77044, 135914,
                                                                       35256, 35346, 80344,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 140939, 0, 3,
                                                                       135914, 77144, 136064,
                                                                       35346, 35436, 80494,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 141164, 0, 3,
                                                                       136064, 77244, 136214,
                                                                       35436, 35526, 80644,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 141389, 0, 3,
                                                                       136214, 77344, 136364,
                                                                       35526, 35616, 80794,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 141614, 0, 3,
                                                                       136364, 77444, 136514,
                                                                       35616, 35706, 80944,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 141839, 0, 3,
                                                                       136514, 77544, 136664,
                                                                       35706, 35796, 81094,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 142064, 0, 3,
                                                                       136664, 77644, 136814,
                                                                       35796, 35886, 81244,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 142289, 0, 3,
                                                                       136814, 77744, 136964,
                                                                       35886, 35976, 81394,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 142514, 0, 3,
                                                                       137114, 77944, 137339,
                                                                       36156, 36282, 81544,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 142829, 0, 3,
                                                                       137339, 78094, 137564,
                                                                       36282, 36408, 81754,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 143144, 0, 3,
                                                                       137564, 78244, 137789,
                                                                       36408, 36534, 81964,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 143459, 0, 3,
                                                                       137789, 78394, 138014,
                                                                       36534, 36660, 82174,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 143774, 0, 3,
                                                                       138014, 78544, 138239,
                                                                       36660, 36786, 82384,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 144089, 0, 3,
                                                                       138239, 78694, 138464,
                                                                       36786, 36912, 82594,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 144404, 0, 3,
                                                                       138464, 78844, 138689,
                                                                       36912, 37038, 82804,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 144719, 0, 3,
                                                                       138689, 78994, 138914,
                                                                       37038, 37164, 83014,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 145034, 0, 3,
                                                                       138914, 79144, 139139,
                                                                       37164, 37290, 83224,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 145349, 0, 3,
                                                                       139139, 79294, 139364,
                                                                       37290, 37416, 83434,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 145664, 0, 3,
                                                                       139364, 79444, 139589,
                                                                       37416, 37542, 83644,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 145979, 0, 3,
                                                                       139814, 79744, 140039,
                                                                       37794, 37920, 83854,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 146294, 0, 3,
                                                                       140039, 79894, 140264,
                                                                       37920, 38046, 84064,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 146609, 0, 3,
                                                                       140264, 80044, 140489,
                                                                       38046, 38172, 84274,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 146924, 0, 3,
                                                                       140489, 80194, 140714,
                                                                       38172, 38298, 84484,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 147239, 0, 3,
                                                                       140714, 80344, 140939,
                                                                       38298, 38424, 84694,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 147554, 0, 3,
                                                                       140939, 80494, 141164,
                                                                       38424, 38550, 84904,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 147869, 0, 3,
                                                                       141164, 80644, 141389,
                                                                       38550, 38676, 85114,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 148184, 0, 3,
                                                                       141389, 80794, 141614,
                                                                       38676, 38802, 85324,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 148499, 0, 3,
                                                                       141614, 80944, 141839,
                                                                       38802, 38928, 85534,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 148814, 0, 3,
                                                                       141839, 81094, 142064,
                                                                       38928, 39054, 85744,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 149129, 0, 3,
                                                                       142064, 81244, 142289,
                                                                       39054, 39180, 85954,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 149444, 0, 3,
                                                                       142514, 81544, 142829,
                                                                       39432, 39600, 86164,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 149864, 0, 3,
                                                                       142829, 81754, 143144,
                                                                       39600, 39768, 86444,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 150284, 0, 3,
                                                                       143144, 81964, 143459,
                                                                       39768, 39936, 86724,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 150704, 0, 3,
                                                                       143459, 82174, 143774,
                                                                       39936, 40104, 87004,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 151124, 0, 3,
                                                                       143774, 82384, 144089,
                                                                       40104, 40272, 87284,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 151544, 0, 3,
                                                                       144089, 82594, 144404,
                                                                       40272, 40440, 87564,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 151964, 0, 3,
                                                                       144404, 82804, 144719,
                                                                       40440, 40608, 87844,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 152384, 0, 3,
                                                                       144719, 83014, 145034,
                                                                       40608, 40776, 88124,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 152804, 0, 3,
                                                                       145034, 83224, 145349,
                                                                       40776, 40944, 88404,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 153224, 0, 3,
                                                                       145349, 83434, 145664,
                                                                       40944, 41112, 88684,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 153644, 0, 3,
                                                                       145979, 83854, 146294,
                                                                       41448, 41616, 88964,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 154064, 0, 3,
                                                                       146294, 84064, 146609,
                                                                       41616, 41784, 89244,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 154484, 0, 3,
                                                                       146609, 84274, 146924,
                                                                       41784, 41952, 89524,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 154904, 0, 3,
                                                                       146924, 84484, 147239,
                                                                       41952, 42120, 89804,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 155324, 0, 3,
                                                                       147239, 84694, 147554,
                                                                       42120, 42288, 90084,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 155744, 0, 3,
                                                                       147554, 84904, 147869,
                                                                       42288, 42456, 90364,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 156164, 0, 3,
                                                                       147869, 85114, 148184,
                                                                       42456, 42624, 90644,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 156584, 0, 3,
                                                                       148184, 85324, 148499,
                                                                       42624, 42792, 90924,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 157004, 0, 3,
                                                                       148499, 85534, 148814,
                                                                       42792, 42960, 91204,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 157424, 0, 3,
                                                                       148814, 85744, 149129,
                                                                       42960, 43128, 91484,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 157844, 0, 3,
                                                                       149444, 86164, 149864,
                                                                       43464, 43680, 91764,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 158384, 0, 3,
                                                                       149864, 86444, 150284,
                                                                       43680, 43896, 92124,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 158924, 0, 3,
                                                                       150284, 86724, 150704,
                                                                       43896, 44112, 92484,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 159464, 0, 3,
                                                                       150704, 87004, 151124,
                                                                       44112, 44328, 92844,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 160004, 0, 3,
                                                                       151124, 87284, 151544,
                                                                       44328, 44544, 93204,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 160544, 0, 3,
                                                                       151544, 87564, 151964,
                                                                       44544, 44760, 93564,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 161084, 0, 3,
                                                                       151964, 87844, 152384,
                                                                       44760, 44976, 93924,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 161624, 0, 3,
                                                                       152384, 88124, 152804,
                                                                       44976, 45192, 94284,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 162164, 0, 3,
                                                                       152804, 88404, 153224,
                                                                       45192, 45408, 94644,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 162704, 0, 3,
                                                                       153644, 88964, 154064,
                                                                       45840, 46056, 95004,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 163244, 0, 3,
                                                                       154064, 89244, 154484,
                                                                       46056, 46272, 95364,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 163784, 0, 3,
                                                                       154484, 89524, 154904,
                                                                       46272, 46488, 95724,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 164324, 0, 3,
                                                                       154904, 89804, 155324,
                                                                       46488, 46704, 96084,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 164864, 0, 3,
                                                                       155324, 90084, 155744,
                                                                       46704, 46920, 96444,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 165404, 0, 3,
                                                                       155744, 90364, 156164,
                                                                       46920, 47136, 96804,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 165944, 0, 3,
                                                                       156164, 90644, 156584,
                                                                       47136, 47352, 97164,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 166484, 0, 3,
                                                                       156584, 90924, 157004,
                                                                       47352, 47568, 97524,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 167024, 0, 3,
                                                                       157004, 91204, 157424,
                                                                       47568, 47784, 97884,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 167564, 0, 3,
                                                                       157844, 91764, 158384,
                                                                       48216, 48486, 98244,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 168239, 0, 3,
                                                                       158384, 92124, 158924,
                                                                       48486, 48756, 98694,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 168914, 0, 3,
                                                                       158924, 92484, 159464,
                                                                       48756, 49026, 99144,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 169589, 0, 3,
                                                                       159464, 92844, 160004,
                                                                       49026, 49296, 99594,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 170264, 0, 3,
                                                                       160004, 93204, 160544,
                                                                       49296, 49566, 100044,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 170939, 0, 3,
                                                                       160544, 93564, 161084,
                                                                       49566, 49836, 100494,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 171614, 0, 3,
                                                                       161084, 93924, 161624,
                                                                       49836, 50106, 100944,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 172289, 0, 3,
                                                                       161624, 94284, 162164,
                                                                       50106, 50376, 101394,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 172964, 0, 3,
                                                                       162704, 95004, 163244,
                                                                       50916, 51186, 101844,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 173639, 0, 3,
                                                                       163244, 95364, 163784,
                                                                       51186, 51456, 102294,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 174314, 0, 3,
                                                                       163784, 95724, 164324,
                                                                       51456, 51726, 102744,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 174989, 0, 3,
                                                                       164324, 96084, 164864,
                                                                       51726, 51996, 103194,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 175664, 0, 3,
                                                                       164864, 96444, 165404,
                                                                       51996, 52266, 103644,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 176339, 0, 3,
                                                                       165404, 96804, 165944,
                                                                       52266, 52536, 104094,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 177014, 0, 3,
                                                                       165944, 97164, 166484,
                                                                       52536, 52806, 104544,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 177689, 0, 3,
                                                                       166484, 97524, 167024,
                                                                       52806, 53076, 104994,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 178364, 0, 3,
                                                                       167564, 98244, 168239,
                                                                       53616, 53946, 105444,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 179189, 0, 3,
                                                                       168239, 98694, 168914,
                                                                       53946, 54276, 105994,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 180014, 0, 3,
                                                                       168914, 99144, 169589,
                                                                       54276, 54606, 106544,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 180839, 0, 3,
                                                                       169589, 99594, 170264,
                                                                       54606, 54936, 107094,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 181664, 0, 3,
                                                                       170264, 100044, 170939,
                                                                       54936, 55266, 107644,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 182489, 0, 3,
                                                                       170939, 100494, 171614,
                                                                       55266, 55596, 108194,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 183314, 0, 3,
                                                                       171614, 100944, 172289,
                                                                       55596, 55926, 108744,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 184139, 0, 3,
                                                                       172964, 101844, 173639,
                                                                       56586, 56916, 109294,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 184964, 0, 3,
                                                                       173639, 102294, 174314,
                                                                       56916, 57246, 109844,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 185789, 0, 3,
                                                                       174314, 102744, 174989,
                                                                       57246, 57576, 110394,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 186614, 0, 3,
                                                                       174989, 103194, 175664,
                                                                       57576, 57906, 110944,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 187439, 0, 3,
                                                                       175664, 103644, 176339,
                                                                       57906, 58236, 111494,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 188264, 0, 3,
                                                                       176339, 104094, 177014,
                                                                       58236, 58566, 112044,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 189089, 0, 3,
                                                                       177014, 104544, 177689,
                                                                       58566, 58896, 112594,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 189914, 0, 3,
                                                                       178364, 105444, 179189,
                                                                       59556, 59952, 113144,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 190904, 0, 3,
                                                                       179189, 105994, 180014,
                                                                       59952, 60348, 113804,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 191894, 0, 3,
                                                                       180014, 106544, 180839,
                                                                       60348, 60744, 114464,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 192884, 0, 3,
                                                                       180839, 107094, 181664,
                                                                       60744, 61140, 115124,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 193874, 0, 3,
                                                                       181664, 107644, 182489,
                                                                       61140, 61536, 115784,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 194864, 0, 3,
                                                                       182489, 108194, 183314,
                                                                       61536, 61932, 116444,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 195854, 0, 3,
                                                                       184139, 109294, 184964,
                                                                       62724, 63120, 117104,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 196844, 0, 3,
                                                                       184964, 109844, 185789,
                                                                       63120, 63516, 117764,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 197834, 0, 3,
                                                                       185789, 110394, 186614,
                                                                       63516, 63912, 118424,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 198824, 0, 3,
                                                                       186614, 110944, 187439,
                                                                       63912, 64308, 119084,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 199814, 0, 3,
                                                                       187439, 111494, 188264,
                                                                       64308, 64704, 119744,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 200804, 0, 3,
                                                                       188264, 112044, 189089,
                                                                       64704, 65100, 120404,
                                                                       ncols, gamma, p, q);

                    compute_prim_sog_three_center_electron_repulsion_0(buffer, 201794, 0, 3,
                                                                       189914, 113144, 190904,
                                                                       65892, 66360, 121064,
                                                                       ncols, gamma, p, q);

                    compute_prim_sog_three_center_electron_repulsion_0(buffer, 202964, 0, 3,
                                                                       190904, 113804, 191894,
                                                                       66360, 66828, 121844,
                                                                       ncols, gamma, p, q);

                    compute_prim_sog_three_center_electron_repulsion_0(buffer, 204134, 0, 3,
                                                                       191894, 114464, 192884,
                                                                       66828, 67296, 122624,
                                                                       ncols, gamma, p, q);

                    compute_prim_sog_three_center_electron_repulsion_0(buffer, 205304, 0, 3,
                                                                       192884, 115124, 193874,
                                                                       67296, 67764, 123404,
                                                                       ncols, gamma, p, q);

                    compute_prim_sog_three_center_electron_repulsion_0(buffer, 206474, 0, 3,
                                                                       193874, 115784, 194864,
                                                                       67764, 68232, 124184,
                                                                       ncols, gamma, p, q);

                    compute_prim_sog_three_center_electron_repulsion_0(buffer, 207644, 0, 3,
                                                                       195854, 117104, 196844,
                                                                       69168, 69636, 124964,
                                                                       ncols, gamma, p, q);

                    compute_prim_sog_three_center_electron_repulsion_0(buffer, 208814, 0, 3,
                                                                       196844, 117764, 197834,
                                                                       69636, 70104, 125744,
                                                                       ncols, gamma, p, q);

                    compute_prim_sog_three_center_electron_repulsion_0(buffer, 209984, 0, 3,
                                                                       197834, 118424, 198824,
                                                                       70104, 70572, 126524,
                                                                       ncols, gamma, p, q);

                    compute_prim_sog_three_center_electron_repulsion_0(buffer, 211154, 0, 3,
                                                                       198824, 119084, 199814,
                                                                       70572, 71040, 127304,
                                                                       ncols, gamma, p, q);

                    compute_prim_sog_three_center_electron_repulsion_0(buffer, 212324, 0, 3,
                                                                       199814, 119744, 200804,
                                                                       71040, 71508, 128084,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 213494, 3, 72444,
                                                                       72454, 128894, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 213515, 3, 72454,
                                                                       72464, 128909, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 213536, 3, 72464,
                                                                       72474, 128924, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 213557, 3, 72474,
                                                                       72484, 128939, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 213578, 3, 72484,
                                                                       72494, 128954, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 213599, 3, 72494,
                                                                       72504, 128969, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 213620, 3, 72504,
                                                                       72514, 128984, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 213641, 3, 72514,
                                                                       72524, 128999, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 213662, 3, 72524,
                                                                       72534, 129014, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 213683, 3, 72534,
                                                                       72544, 129029, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 213704, 3, 72544,
                                                                       72554, 129044, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 213725, 3, 72554,
                                                                       72564, 129059, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 213746, 3, 72564,
                                                                       72574, 129074, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 213767, 3, 72574,
                                                                       72584, 129089, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 213788, 3, 72604,
                                                                       72614, 129134, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 213809, 3, 72614,
                                                                       72624, 129149, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 213830, 3, 72624,
                                                                       72634, 129164, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 213851, 3, 72634,
                                                                       72644, 129179, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 213872, 3, 72644,
                                                                       72654, 129194, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 213893, 3, 72654,
                                                                       72664, 129209, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 213914, 3, 72664,
                                                                       72674, 129224, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 213935, 3, 72674,
                                                                       72684, 129239, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 213956, 3, 72684,
                                                                       72694, 129254, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 213977, 3, 72694,
                                                                       72704, 129269, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 213998, 3, 72704,
                                                                       72714, 129284, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 214019, 3, 72714,
                                                                       72724, 129299, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 214040, 3, 72724,
                                                                       72734, 129314, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 214061, 3, 72734,
                                                                       72744, 129329, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 214082, 0, 3,
                                                                       213494, 128894, 213515,
                                                                       72764, 72794, 129434,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 214145, 0, 3,
                                                                       213515, 128909, 213536,
                                                                       72794, 72824, 129479,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 214208, 0, 3,
                                                                       213536, 128924, 213557,
                                                                       72824, 72854, 129524,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 214271, 0, 3,
                                                                       213557, 128939, 213578,
                                                                       72854, 72884, 129569,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 214334, 0, 3,
                                                                       213578, 128954, 213599,
                                                                       72884, 72914, 129614,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 214397, 0, 3,
                                                                       213599, 128969, 213620,
                                                                       72914, 72944, 129659,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 214460, 0, 3,
                                                                       213620, 128984, 213641,
                                                                       72944, 72974, 129704,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 214523, 0, 3,
                                                                       213641, 128999, 213662,
                                                                       72974, 73004, 129749,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 214586, 0, 3,
                                                                       213662, 129014, 213683,
                                                                       73004, 73034, 129794,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 214649, 0, 3,
                                                                       213683, 129029, 213704,
                                                                       73034, 73064, 129839,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 214712, 0, 3,
                                                                       213704, 129044, 213725,
                                                                       73064, 73094, 129884,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 214775, 0, 3,
                                                                       213725, 129059, 213746,
                                                                       73094, 73124, 129929,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 214838, 0, 3,
                                                                       213746, 129074, 213767,
                                                                       73124, 73154, 129974,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 214901, 0, 3,
                                                                       213788, 129134, 213809,
                                                                       73214, 73244, 130109,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 214964, 0, 3,
                                                                       213809, 129149, 213830,
                                                                       73244, 73274, 130154,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 215027, 0, 3,
                                                                       213830, 129164, 213851,
                                                                       73274, 73304, 130199,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 215090, 0, 3,
                                                                       213851, 129179, 213872,
                                                                       73304, 73334, 130244,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 215153, 0, 3,
                                                                       213872, 129194, 213893,
                                                                       73334, 73364, 130289,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 215216, 0, 3,
                                                                       213893, 129209, 213914,
                                                                       73364, 73394, 130334,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 215279, 0, 3,
                                                                       213914, 129224, 213935,
                                                                       73394, 73424, 130379,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 215342, 0, 3,
                                                                       213935, 129239, 213956,
                                                                       73424, 73454, 130424,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 215405, 0, 3,
                                                                       213956, 129254, 213977,
                                                                       73454, 73484, 130469,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 215468, 0, 3,
                                                                       213977, 129269, 213998,
                                                                       73484, 73514, 130514,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 215531, 0, 3,
                                                                       213998, 129284, 214019,
                                                                       73514, 73544, 130559,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 215594, 0, 3,
                                                                       214019, 129299, 214040,
                                                                       73544, 73574, 130604,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 215657, 0, 3,
                                                                       214040, 129314, 214061,
                                                                       73574, 73604, 130649,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 215720, 0, 3,
                                                                       214082, 129434, 214145,
                                                                       73664, 73724, 130874,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 215846, 0, 3,
                                                                       214145, 129479, 214208,
                                                                       73724, 73784, 130964,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 215972, 0, 3,
                                                                       214208, 129524, 214271,
                                                                       73784, 73844, 131054,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 216098, 0, 3,
                                                                       214271, 129569, 214334,
                                                                       73844, 73904, 131144,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 216224, 0, 3,
                                                                       214334, 129614, 214397,
                                                                       73904, 73964, 131234,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 216350, 0, 3,
                                                                       214397, 129659, 214460,
                                                                       73964, 74024, 131324,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 216476, 0, 3,
                                                                       214460, 129704, 214523,
                                                                       74024, 74084, 131414,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 216602, 0, 3,
                                                                       214523, 129749, 214586,
                                                                       74084, 74144, 131504,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 216728, 0, 3,
                                                                       214586, 129794, 214649,
                                                                       74144, 74204, 131594,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 216854, 0, 3,
                                                                       214649, 129839, 214712,
                                                                       74204, 74264, 131684,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 216980, 0, 3,
                                                                       214712, 129884, 214775,
                                                                       74264, 74324, 131774,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 217106, 0, 3,
                                                                       214775, 129929, 214838,
                                                                       74324, 74384, 131864,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 217232, 0, 3,
                                                                       214901, 130109, 214964,
                                                                       74504, 74564, 132134,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 217358, 0, 3,
                                                                       214964, 130154, 215027,
                                                                       74564, 74624, 132224,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 217484, 0, 3,
                                                                       215027, 130199, 215090,
                                                                       74624, 74684, 132314,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 217610, 0, 3,
                                                                       215090, 130244, 215153,
                                                                       74684, 74744, 132404,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 217736, 0, 3,
                                                                       215153, 130289, 215216,
                                                                       74744, 74804, 132494,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 217862, 0, 3,
                                                                       215216, 130334, 215279,
                                                                       74804, 74864, 132584,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 217988, 0, 3,
                                                                       215279, 130379, 215342,
                                                                       74864, 74924, 132674,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 218114, 0, 3,
                                                                       215342, 130424, 215405,
                                                                       74924, 74984, 132764,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 218240, 0, 3,
                                                                       215405, 130469, 215468,
                                                                       74984, 75044, 132854,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 218366, 0, 3,
                                                                       215468, 130514, 215531,
                                                                       75044, 75104, 132944,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 218492, 0, 3,
                                                                       215531, 130559, 215594,
                                                                       75104, 75164, 133034,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 218618, 0, 3,
                                                                       215594, 130604, 215657,
                                                                       75164, 75224, 133124,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 218744, 0, 3,
                                                                       215720, 130874, 215846,
                                                                       75344, 75444, 133514,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 218954, 0, 3,
                                                                       215846, 130964, 215972,
                                                                       75444, 75544, 133664,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 219164, 0, 3,
                                                                       215972, 131054, 216098,
                                                                       75544, 75644, 133814,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 219374, 0, 3,
                                                                       216098, 131144, 216224,
                                                                       75644, 75744, 133964,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 219584, 0, 3,
                                                                       216224, 131234, 216350,
                                                                       75744, 75844, 134114,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 219794, 0, 3,
                                                                       216350, 131324, 216476,
                                                                       75844, 75944, 134264,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 220004, 0, 3,
                                                                       216476, 131414, 216602,
                                                                       75944, 76044, 134414,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 220214, 0, 3,
                                                                       216602, 131504, 216728,
                                                                       76044, 76144, 134564,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 220424, 0, 3,
                                                                       216728, 131594, 216854,
                                                                       76144, 76244, 134714,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 220634, 0, 3,
                                                                       216854, 131684, 216980,
                                                                       76244, 76344, 134864,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 220844, 0, 3,
                                                                       216980, 131774, 217106,
                                                                       76344, 76444, 135014,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 221054, 0, 3,
                                                                       217232, 132134, 217358,
                                                                       76644, 76744, 135464,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 221264, 0, 3,
                                                                       217358, 132224, 217484,
                                                                       76744, 76844, 135614,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 221474, 0, 3,
                                                                       217484, 132314, 217610,
                                                                       76844, 76944, 135764,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 221684, 0, 3,
                                                                       217610, 132404, 217736,
                                                                       76944, 77044, 135914,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 221894, 0, 3,
                                                                       217736, 132494, 217862,
                                                                       77044, 77144, 136064,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 222104, 0, 3,
                                                                       217862, 132584, 217988,
                                                                       77144, 77244, 136214,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 222314, 0, 3,
                                                                       217988, 132674, 218114,
                                                                       77244, 77344, 136364,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 222524, 0, 3,
                                                                       218114, 132764, 218240,
                                                                       77344, 77444, 136514,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 222734, 0, 3,
                                                                       218240, 132854, 218366,
                                                                       77444, 77544, 136664,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 222944, 0, 3,
                                                                       218366, 132944, 218492,
                                                                       77544, 77644, 136814,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 223154, 0, 3,
                                                                       218492, 133034, 218618,
                                                                       77644, 77744, 136964,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 223364, 0, 3,
                                                                       218744, 133514, 218954,
                                                                       77944, 78094, 137564,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 223679, 0, 3,
                                                                       218954, 133664, 219164,
                                                                       78094, 78244, 137789,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 223994, 0, 3,
                                                                       219164, 133814, 219374,
                                                                       78244, 78394, 138014,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 224309, 0, 3,
                                                                       219374, 133964, 219584,
                                                                       78394, 78544, 138239,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 224624, 0, 3,
                                                                       219584, 134114, 219794,
                                                                       78544, 78694, 138464,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 224939, 0, 3,
                                                                       219794, 134264, 220004,
                                                                       78694, 78844, 138689,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 225254, 0, 3,
                                                                       220004, 134414, 220214,
                                                                       78844, 78994, 138914,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 225569, 0, 3,
                                                                       220214, 134564, 220424,
                                                                       78994, 79144, 139139,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 225884, 0, 3,
                                                                       220424, 134714, 220634,
                                                                       79144, 79294, 139364,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 226199, 0, 3,
                                                                       220634, 134864, 220844,
                                                                       79294, 79444, 139589,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 226514, 0, 3,
                                                                       221054, 135464, 221264,
                                                                       79744, 79894, 140264,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 226829, 0, 3,
                                                                       221264, 135614, 221474,
                                                                       79894, 80044, 140489,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 227144, 0, 3,
                                                                       221474, 135764, 221684,
                                                                       80044, 80194, 140714,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 227459, 0, 3,
                                                                       221684, 135914, 221894,
                                                                       80194, 80344, 140939,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 227774, 0, 3,
                                                                       221894, 136064, 222104,
                                                                       80344, 80494, 141164,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 228089, 0, 3,
                                                                       222104, 136214, 222314,
                                                                       80494, 80644, 141389,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 228404, 0, 3,
                                                                       222314, 136364, 222524,
                                                                       80644, 80794, 141614,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 228719, 0, 3,
                                                                       222524, 136514, 222734,
                                                                       80794, 80944, 141839,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 229034, 0, 3,
                                                                       222734, 136664, 222944,
                                                                       80944, 81094, 142064,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 229349, 0, 3,
                                                                       222944, 136814, 223154,
                                                                       81094, 81244, 142289,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 229664, 0, 3,
                                                                       223364, 137564, 223679,
                                                                       81544, 81754, 143144,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 230105, 0, 3,
                                                                       223679, 137789, 223994,
                                                                       81754, 81964, 143459,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 230546, 0, 3,
                                                                       223994, 138014, 224309,
                                                                       81964, 82174, 143774,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 230987, 0, 3,
                                                                       224309, 138239, 224624,
                                                                       82174, 82384, 144089,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 231428, 0, 3,
                                                                       224624, 138464, 224939,
                                                                       82384, 82594, 144404,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 231869, 0, 3,
                                                                       224939, 138689, 225254,
                                                                       82594, 82804, 144719,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 232310, 0, 3,
                                                                       225254, 138914, 225569,
                                                                       82804, 83014, 145034,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 232751, 0, 3,
                                                                       225569, 139139, 225884,
                                                                       83014, 83224, 145349,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 233192, 0, 3,
                                                                       225884, 139364, 226199,
                                                                       83224, 83434, 145664,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 233633, 0, 3,
                                                                       226514, 140264, 226829,
                                                                       83854, 84064, 146609,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 234074, 0, 3,
                                                                       226829, 140489, 227144,
                                                                       84064, 84274, 146924,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 234515, 0, 3,
                                                                       227144, 140714, 227459,
                                                                       84274, 84484, 147239,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 234956, 0, 3,
                                                                       227459, 140939, 227774,
                                                                       84484, 84694, 147554,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 235397, 0, 3,
                                                                       227774, 141164, 228089,
                                                                       84694, 84904, 147869,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 235838, 0, 3,
                                                                       228089, 141389, 228404,
                                                                       84904, 85114, 148184,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 236279, 0, 3,
                                                                       228404, 141614, 228719,
                                                                       85114, 85324, 148499,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 236720, 0, 3,
                                                                       228719, 141839, 229034,
                                                                       85324, 85534, 148814,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 237161, 0, 3,
                                                                       229034, 142064, 229349,
                                                                       85534, 85744, 149129,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 237602, 0, 3,
                                                                       229664, 143144, 230105,
                                                                       86164, 86444, 150284,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 238190, 0, 3,
                                                                       230105, 143459, 230546,
                                                                       86444, 86724, 150704,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 238778, 0, 3,
                                                                       230546, 143774, 230987,
                                                                       86724, 87004, 151124,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 239366, 0, 3,
                                                                       230987, 144089, 231428,
                                                                       87004, 87284, 151544,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 239954, 0, 3,
                                                                       231428, 144404, 231869,
                                                                       87284, 87564, 151964,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 240542, 0, 3,
                                                                       231869, 144719, 232310,
                                                                       87564, 87844, 152384,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 241130, 0, 3,
                                                                       232310, 145034, 232751,
                                                                       87844, 88124, 152804,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 241718, 0, 3,
                                                                       232751, 145349, 233192,
                                                                       88124, 88404, 153224,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 242306, 0, 3,
                                                                       233633, 146609, 234074,
                                                                       88964, 89244, 154484,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 242894, 0, 3,
                                                                       234074, 146924, 234515,
                                                                       89244, 89524, 154904,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 243482, 0, 3,
                                                                       234515, 147239, 234956,
                                                                       89524, 89804, 155324,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 244070, 0, 3,
                                                                       234956, 147554, 235397,
                                                                       89804, 90084, 155744,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 244658, 0, 3,
                                                                       235397, 147869, 235838,
                                                                       90084, 90364, 156164,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 245246, 0, 3,
                                                                       235838, 148184, 236279,
                                                                       90364, 90644, 156584,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 245834, 0, 3,
                                                                       236279, 148499, 236720,
                                                                       90644, 90924, 157004,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 246422, 0, 3,
                                                                       236720, 148814, 237161,
                                                                       90924, 91204, 157424,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 247010, 0, 3,
                                                                       237602, 150284, 238190,
                                                                       91764, 92124, 158924,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 247766, 0, 3,
                                                                       238190, 150704, 238778,
                                                                       92124, 92484, 159464,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 248522, 0, 3,
                                                                       238778, 151124, 239366,
                                                                       92484, 92844, 160004,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 249278, 0, 3,
                                                                       239366, 151544, 239954,
                                                                       92844, 93204, 160544,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 250034, 0, 3,
                                                                       239954, 151964, 240542,
                                                                       93204, 93564, 161084,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 250790, 0, 3,
                                                                       240542, 152384, 241130,
                                                                       93564, 93924, 161624,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 251546, 0, 3,
                                                                       241130, 152804, 241718,
                                                                       93924, 94284, 162164,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 252302, 0, 3,
                                                                       242306, 154484, 242894,
                                                                       95004, 95364, 163784,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 253058, 0, 3,
                                                                       242894, 154904, 243482,
                                                                       95364, 95724, 164324,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 253814, 0, 3,
                                                                       243482, 155324, 244070,
                                                                       95724, 96084, 164864,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 254570, 0, 3,
                                                                       244070, 155744, 244658,
                                                                       96084, 96444, 165404,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 255326, 0, 3,
                                                                       244658, 156164, 245246,
                                                                       96444, 96804, 165944,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 256082, 0, 3,
                                                                       245246, 156584, 245834,
                                                                       96804, 97164, 166484,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 256838, 0, 3,
                                                                       245834, 157004, 246422,
                                                                       97164, 97524, 167024,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 257594, 0, 3,
                                                                       247010, 158924, 247766,
                                                                       98244, 98694, 168914,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 258539, 0, 3,
                                                                       247766, 159464, 248522,
                                                                       98694, 99144, 169589,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 259484, 0, 3,
                                                                       248522, 160004, 249278,
                                                                       99144, 99594, 170264,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 260429, 0, 3,
                                                                       249278, 160544, 250034,
                                                                       99594, 100044, 170939,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 261374, 0, 3,
                                                                       250034, 161084, 250790,
                                                                       100044, 100494, 171614,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 262319, 0, 3,
                                                                       250790, 161624, 251546,
                                                                       100494, 100944, 172289,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 263264, 0, 3,
                                                                       252302, 163784, 253058,
                                                                       101844, 102294, 174314,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 264209, 0, 3,
                                                                       253058, 164324, 253814,
                                                                       102294, 102744, 174989,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 265154, 0, 3,
                                                                       253814, 164864, 254570,
                                                                       102744, 103194, 175664,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 266099, 0, 3,
                                                                       254570, 165404, 255326,
                                                                       103194, 103644, 176339,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 267044, 0, 3,
                                                                       255326, 165944, 256082,
                                                                       103644, 104094, 177014,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 267989, 0, 3,
                                                                       256082, 166484, 256838,
                                                                       104094, 104544, 177689,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 268934, 0, 3,
                                                                       257594, 168914, 258539,
                                                                       105444, 105994, 180014,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 270089, 0, 3,
                                                                       258539, 169589, 259484,
                                                                       105994, 106544, 180839,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 271244, 0, 3,
                                                                       259484, 170264, 260429,
                                                                       106544, 107094, 181664,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 272399, 0, 3,
                                                                       260429, 170939, 261374,
                                                                       107094, 107644, 182489,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 273554, 0, 3,
                                                                       261374, 171614, 262319,
                                                                       107644, 108194, 183314,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 274709, 0, 3,
                                                                       263264, 174314, 264209,
                                                                       109294, 109844, 185789,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 275864, 0, 3,
                                                                       264209, 174989, 265154,
                                                                       109844, 110394, 186614,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 277019, 0, 3,
                                                                       265154, 175664, 266099,
                                                                       110394, 110944, 187439,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 278174, 0, 3,
                                                                       266099, 176339, 267044,
                                                                       110944, 111494, 188264,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 279329, 0, 3,
                                                                       267044, 177014, 267989,
                                                                       111494, 112044, 189089,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 280484, 0, 3,
                                                                       268934, 180014, 270089,
                                                                       113144, 113804, 191894,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 281870, 0, 3,
                                                                       270089, 180839, 271244,
                                                                       113804, 114464, 192884,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 283256, 0, 3,
                                                                       271244, 181664, 272399,
                                                                       114464, 115124, 193874,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 284642, 0, 3,
                                                                       272399, 182489, 273554,
                                                                       115124, 115784, 194864,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 286028, 0, 3,
                                                                       274709, 185789, 275864,
                                                                       117104, 117764, 197834,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 287414, 0, 3,
                                                                       275864, 186614, 277019,
                                                                       117764, 118424, 198824,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 288800, 0, 3,
                                                                       277019, 187439, 278174,
                                                                       118424, 119084, 199814,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 290186, 0, 3,
                                                                       278174, 188264, 279329,
                                                                       119084, 119744, 200804,
                                                                       ncols, gamma, p, q);

                    compute_prim_soh_three_center_electron_repulsion_0(buffer, 291572, 0, 3,
                                                                       280484, 191894, 281870,
                                                                       121064, 121844, 204134,
                                                                       ncols, gamma, p, q);

                    compute_prim_soh_three_center_electron_repulsion_0(buffer, 293210, 0, 3,
                                                                       281870, 192884, 283256,
                                                                       121844, 122624, 205304,
                                                                       ncols, gamma, p, q);

                    compute_prim_soh_three_center_electron_repulsion_0(buffer, 294848, 0, 3,
                                                                       283256, 193874, 284642,
                                                                       122624, 123404, 206474,
                                                                       ncols, gamma, p, q);

                    compute_prim_soh_three_center_electron_repulsion_0(buffer, 296486, 0, 3,
                                                                       286028, 197834, 287414,
                                                                       124964, 125744, 209984,
                                                                       ncols, gamma, p, q);

                    compute_prim_soh_three_center_electron_repulsion_0(buffer, 298124, 0, 3,
                                                                       287414, 198824, 288800,
                                                                       125744, 126524, 211154,
                                                                       ncols, gamma, p, q);

                    compute_prim_soh_three_center_electron_repulsion_0(buffer, 299762, 0, 3,
                                                                       288800, 199814, 290186,
                                                                       126524, 127304, 212324,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 301400, 3, 128864,
                                                                       128879, 213494, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 301428, 3, 128879,
                                                                       128894, 213515, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 301456, 3, 128894,
                                                                       128909, 213536, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 301484, 3, 128909,
                                                                       128924, 213557, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 301512, 3, 128924,
                                                                       128939, 213578, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 301540, 3, 128939,
                                                                       128954, 213599, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 301568, 3, 128954,
                                                                       128969, 213620, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 301596, 3, 128969,
                                                                       128984, 213641, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 301624, 3, 128984,
                                                                       128999, 213662, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 301652, 3, 128999,
                                                                       129014, 213683, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 301680, 3, 129014,
                                                                       129029, 213704, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 301708, 3, 129029,
                                                                       129044, 213725, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 301736, 3, 129044,
                                                                       129059, 213746, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 301764, 3, 129059,
                                                                       129074, 213767, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 301792, 3, 129104,
                                                                       129119, 213788, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 301820, 3, 129119,
                                                                       129134, 213809, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 301848, 3, 129134,
                                                                       129149, 213830, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 301876, 3, 129149,
                                                                       129164, 213851, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 301904, 3, 129164,
                                                                       129179, 213872, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 301932, 3, 129179,
                                                                       129194, 213893, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 301960, 3, 129194,
                                                                       129209, 213914, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 301988, 3, 129209,
                                                                       129224, 213935, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 302016, 3, 129224,
                                                                       129239, 213956, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 302044, 3, 129239,
                                                                       129254, 213977, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 302072, 3, 129254,
                                                                       129269, 213998, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 302100, 3, 129269,
                                                                       129284, 214019, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 302128, 3, 129284,
                                                                       129299, 214040, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 302156, 3, 129299,
                                                                       129314, 214061, ncols,
                                                                       gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 302184, 0, 3,
                                                                       301400, 213494, 301428,
                                                                       129344, 129389, 214082,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 302268, 0, 3,
                                                                       301428, 213515, 301456,
                                                                       129389, 129434, 214145,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 302352, 0, 3,
                                                                       301456, 213536, 301484,
                                                                       129434, 129479, 214208,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 302436, 0, 3,
                                                                       301484, 213557, 301512,
                                                                       129479, 129524, 214271,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 302520, 0, 3,
                                                                       301512, 213578, 301540,
                                                                       129524, 129569, 214334,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 302604, 0, 3,
                                                                       301540, 213599, 301568,
                                                                       129569, 129614, 214397,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 302688, 0, 3,
                                                                       301568, 213620, 301596,
                                                                       129614, 129659, 214460,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 302772, 0, 3,
                                                                       301596, 213641, 301624,
                                                                       129659, 129704, 214523,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 302856, 0, 3,
                                                                       301624, 213662, 301652,
                                                                       129704, 129749, 214586,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 302940, 0, 3,
                                                                       301652, 213683, 301680,
                                                                       129749, 129794, 214649,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 303024, 0, 3,
                                                                       301680, 213704, 301708,
                                                                       129794, 129839, 214712,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 303108, 0, 3,
                                                                       301708, 213725, 301736,
                                                                       129839, 129884, 214775,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 303192, 0, 3,
                                                                       301736, 213746, 301764,
                                                                       129884, 129929, 214838,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 303276, 0, 3,
                                                                       301792, 213788, 301820,
                                                                       130019, 130064, 214901,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 303360, 0, 3,
                                                                       301820, 213809, 301848,
                                                                       130064, 130109, 214964,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 303444, 0, 3,
                                                                       301848, 213830, 301876,
                                                                       130109, 130154, 215027,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 303528, 0, 3,
                                                                       301876, 213851, 301904,
                                                                       130154, 130199, 215090,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 303612, 0, 3,
                                                                       301904, 213872, 301932,
                                                                       130199, 130244, 215153,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 303696, 0, 3,
                                                                       301932, 213893, 301960,
                                                                       130244, 130289, 215216,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 303780, 0, 3,
                                                                       301960, 213914, 301988,
                                                                       130289, 130334, 215279,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 303864, 0, 3,
                                                                       301988, 213935, 302016,
                                                                       130334, 130379, 215342,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 303948, 0, 3,
                                                                       302016, 213956, 302044,
                                                                       130379, 130424, 215405,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 304032, 0, 3,
                                                                       302044, 213977, 302072,
                                                                       130424, 130469, 215468,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 304116, 0, 3,
                                                                       302072, 213998, 302100,
                                                                       130469, 130514, 215531,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 304200, 0, 3,
                                                                       302100, 214019, 302128,
                                                                       130514, 130559, 215594,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 304284, 0, 3,
                                                                       302128, 214040, 302156,
                                                                       130559, 130604, 215657,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 304368, 0, 3,
                                                                       302184, 214082, 302268,
                                                                       130694, 130784, 215720,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 304536, 0, 3,
                                                                       302268, 214145, 302352,
                                                                       130784, 130874, 215846,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 304704, 0, 3,
                                                                       302352, 214208, 302436,
                                                                       130874, 130964, 215972,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 304872, 0, 3,
                                                                       302436, 214271, 302520,
                                                                       130964, 131054, 216098,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 305040, 0, 3,
                                                                       302520, 214334, 302604,
                                                                       131054, 131144, 216224,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 305208, 0, 3,
                                                                       302604, 214397, 302688,
                                                                       131144, 131234, 216350,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 305376, 0, 3,
                                                                       302688, 214460, 302772,
                                                                       131234, 131324, 216476,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 305544, 0, 3,
                                                                       302772, 214523, 302856,
                                                                       131324, 131414, 216602,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 305712, 0, 3,
                                                                       302856, 214586, 302940,
                                                                       131414, 131504, 216728,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 305880, 0, 3,
                                                                       302940, 214649, 303024,
                                                                       131504, 131594, 216854,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 306048, 0, 3,
                                                                       303024, 214712, 303108,
                                                                       131594, 131684, 216980,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 306216, 0, 3,
                                                                       303108, 214775, 303192,
                                                                       131684, 131774, 217106,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 306384, 0, 3,
                                                                       303276, 214901, 303360,
                                                                       131954, 132044, 217232,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 306552, 0, 3,
                                                                       303360, 214964, 303444,
                                                                       132044, 132134, 217358,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 306720, 0, 3,
                                                                       303444, 215027, 303528,
                                                                       132134, 132224, 217484,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 306888, 0, 3,
                                                                       303528, 215090, 303612,
                                                                       132224, 132314, 217610,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 307056, 0, 3,
                                                                       303612, 215153, 303696,
                                                                       132314, 132404, 217736,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 307224, 0, 3,
                                                                       303696, 215216, 303780,
                                                                       132404, 132494, 217862,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 307392, 0, 3,
                                                                       303780, 215279, 303864,
                                                                       132494, 132584, 217988,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 307560, 0, 3,
                                                                       303864, 215342, 303948,
                                                                       132584, 132674, 218114,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 307728, 0, 3,
                                                                       303948, 215405, 304032,
                                                                       132674, 132764, 218240,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 307896, 0, 3,
                                                                       304032, 215468, 304116,
                                                                       132764, 132854, 218366,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 308064, 0, 3,
                                                                       304116, 215531, 304200,
                                                                       132854, 132944, 218492,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 308232, 0, 3,
                                                                       304200, 215594, 304284,
                                                                       132944, 133034, 218618,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 308400, 0, 3,
                                                                       304368, 215720, 304536,
                                                                       133214, 133364, 218744,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 308680, 0, 3,
                                                                       304536, 215846, 304704,
                                                                       133364, 133514, 218954,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 308960, 0, 3,
                                                                       304704, 215972, 304872,
                                                                       133514, 133664, 219164,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 309240, 0, 3,
                                                                       304872, 216098, 305040,
                                                                       133664, 133814, 219374,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 309520, 0, 3,
                                                                       305040, 216224, 305208,
                                                                       133814, 133964, 219584,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 309800, 0, 3,
                                                                       305208, 216350, 305376,
                                                                       133964, 134114, 219794,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 310080, 0, 3,
                                                                       305376, 216476, 305544,
                                                                       134114, 134264, 220004,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 310360, 0, 3,
                                                                       305544, 216602, 305712,
                                                                       134264, 134414, 220214,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 310640, 0, 3,
                                                                       305712, 216728, 305880,
                                                                       134414, 134564, 220424,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 310920, 0, 3,
                                                                       305880, 216854, 306048,
                                                                       134564, 134714, 220634,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 311200, 0, 3,
                                                                       306048, 216980, 306216,
                                                                       134714, 134864, 220844,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 311480, 0, 3,
                                                                       306384, 217232, 306552,
                                                                       135164, 135314, 221054,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 311760, 0, 3,
                                                                       306552, 217358, 306720,
                                                                       135314, 135464, 221264,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 312040, 0, 3,
                                                                       306720, 217484, 306888,
                                                                       135464, 135614, 221474,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 312320, 0, 3,
                                                                       306888, 217610, 307056,
                                                                       135614, 135764, 221684,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 312600, 0, 3,
                                                                       307056, 217736, 307224,
                                                                       135764, 135914, 221894,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 312880, 0, 3,
                                                                       307224, 217862, 307392,
                                                                       135914, 136064, 222104,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 313160, 0, 3,
                                                                       307392, 217988, 307560,
                                                                       136064, 136214, 222314,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 313440, 0, 3,
                                                                       307560, 218114, 307728,
                                                                       136214, 136364, 222524,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 313720, 0, 3,
                                                                       307728, 218240, 307896,
                                                                       136364, 136514, 222734,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 314000, 0, 3,
                                                                       307896, 218366, 308064,
                                                                       136514, 136664, 222944,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 314280, 0, 3,
                                                                       308064, 218492, 308232,
                                                                       136664, 136814, 223154,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 314560, 0, 3,
                                                                       308400, 218744, 308680,
                                                                       137114, 137339, 223364,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 314980, 0, 3,
                                                                       308680, 218954, 308960,
                                                                       137339, 137564, 223679,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 315400, 0, 3,
                                                                       308960, 219164, 309240,
                                                                       137564, 137789, 223994,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 315820, 0, 3,
                                                                       309240, 219374, 309520,
                                                                       137789, 138014, 224309,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 316240, 0, 3,
                                                                       309520, 219584, 309800,
                                                                       138014, 138239, 224624,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 316660, 0, 3,
                                                                       309800, 219794, 310080,
                                                                       138239, 138464, 224939,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 317080, 0, 3,
                                                                       310080, 220004, 310360,
                                                                       138464, 138689, 225254,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 317500, 0, 3,
                                                                       310360, 220214, 310640,
                                                                       138689, 138914, 225569,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 317920, 0, 3,
                                                                       310640, 220424, 310920,
                                                                       138914, 139139, 225884,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 318340, 0, 3,
                                                                       310920, 220634, 311200,
                                                                       139139, 139364, 226199,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 318760, 0, 3,
                                                                       311480, 221054, 311760,
                                                                       139814, 140039, 226514,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 319180, 0, 3,
                                                                       311760, 221264, 312040,
                                                                       140039, 140264, 226829,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 319600, 0, 3,
                                                                       312040, 221474, 312320,
                                                                       140264, 140489, 227144,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 320020, 0, 3,
                                                                       312320, 221684, 312600,
                                                                       140489, 140714, 227459,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 320440, 0, 3,
                                                                       312600, 221894, 312880,
                                                                       140714, 140939, 227774,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 320860, 0, 3,
                                                                       312880, 222104, 313160,
                                                                       140939, 141164, 228089,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 321280, 0, 3,
                                                                       313160, 222314, 313440,
                                                                       141164, 141389, 228404,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 321700, 0, 3,
                                                                       313440, 222524, 313720,
                                                                       141389, 141614, 228719,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 322120, 0, 3,
                                                                       313720, 222734, 314000,
                                                                       141614, 141839, 229034,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 322540, 0, 3,
                                                                       314000, 222944, 314280,
                                                                       141839, 142064, 229349,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 322960, 0, 3,
                                                                       314560, 223364, 314980,
                                                                       142514, 142829, 229664,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 323548, 0, 3,
                                                                       314980, 223679, 315400,
                                                                       142829, 143144, 230105,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 324136, 0, 3,
                                                                       315400, 223994, 315820,
                                                                       143144, 143459, 230546,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 324724, 0, 3,
                                                                       315820, 224309, 316240,
                                                                       143459, 143774, 230987,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 325312, 0, 3,
                                                                       316240, 224624, 316660,
                                                                       143774, 144089, 231428,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 325900, 0, 3,
                                                                       316660, 224939, 317080,
                                                                       144089, 144404, 231869,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 326488, 0, 3,
                                                                       317080, 225254, 317500,
                                                                       144404, 144719, 232310,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 327076, 0, 3,
                                                                       317500, 225569, 317920,
                                                                       144719, 145034, 232751,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 327664, 0, 3,
                                                                       317920, 225884, 318340,
                                                                       145034, 145349, 233192,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 328252, 0, 3,
                                                                       318760, 226514, 319180,
                                                                       145979, 146294, 233633,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 328840, 0, 3,
                                                                       319180, 226829, 319600,
                                                                       146294, 146609, 234074,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 329428, 0, 3,
                                                                       319600, 227144, 320020,
                                                                       146609, 146924, 234515,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 330016, 0, 3,
                                                                       320020, 227459, 320440,
                                                                       146924, 147239, 234956,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 330604, 0, 3,
                                                                       320440, 227774, 320860,
                                                                       147239, 147554, 235397,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 331192, 0, 3,
                                                                       320860, 228089, 321280,
                                                                       147554, 147869, 235838,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 331780, 0, 3,
                                                                       321280, 228404, 321700,
                                                                       147869, 148184, 236279,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 332368, 0, 3,
                                                                       321700, 228719, 322120,
                                                                       148184, 148499, 236720,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 332956, 0, 3,
                                                                       322120, 229034, 322540,
                                                                       148499, 148814, 237161,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 333544, 0, 3,
                                                                       322960, 229664, 323548,
                                                                       149444, 149864, 237602,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 334328, 0, 3,
                                                                       323548, 230105, 324136,
                                                                       149864, 150284, 238190,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 335112, 0, 3,
                                                                       324136, 230546, 324724,
                                                                       150284, 150704, 238778,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 335896, 0, 3,
                                                                       324724, 230987, 325312,
                                                                       150704, 151124, 239366,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 336680, 0, 3,
                                                                       325312, 231428, 325900,
                                                                       151124, 151544, 239954,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 337464, 0, 3,
                                                                       325900, 231869, 326488,
                                                                       151544, 151964, 240542,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 338248, 0, 3,
                                                                       326488, 232310, 327076,
                                                                       151964, 152384, 241130,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 339032, 0, 3,
                                                                       327076, 232751, 327664,
                                                                       152384, 152804, 241718,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 339816, 0, 3,
                                                                       328252, 233633, 328840,
                                                                       153644, 154064, 242306,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 340600, 0, 3,
                                                                       328840, 234074, 329428,
                                                                       154064, 154484, 242894,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 341384, 0, 3,
                                                                       329428, 234515, 330016,
                                                                       154484, 154904, 243482,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 342168, 0, 3,
                                                                       330016, 234956, 330604,
                                                                       154904, 155324, 244070,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 342952, 0, 3,
                                                                       330604, 235397, 331192,
                                                                       155324, 155744, 244658,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 343736, 0, 3,
                                                                       331192, 235838, 331780,
                                                                       155744, 156164, 245246,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 344520, 0, 3,
                                                                       331780, 236279, 332368,
                                                                       156164, 156584, 245834,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 345304, 0, 3,
                                                                       332368, 236720, 332956,
                                                                       156584, 157004, 246422,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 346088, 0, 3,
                                                                       333544, 237602, 334328,
                                                                       157844, 158384, 247010,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 347096, 0, 3,
                                                                       334328, 238190, 335112,
                                                                       158384, 158924, 247766,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 348104, 0, 3,
                                                                       335112, 238778, 335896,
                                                                       158924, 159464, 248522,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 349112, 0, 3,
                                                                       335896, 239366, 336680,
                                                                       159464, 160004, 249278,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 350120, 0, 3,
                                                                       336680, 239954, 337464,
                                                                       160004, 160544, 250034,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 351128, 0, 3,
                                                                       337464, 240542, 338248,
                                                                       160544, 161084, 250790,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 352136, 0, 3,
                                                                       338248, 241130, 339032,
                                                                       161084, 161624, 251546,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 353144, 0, 3,
                                                                       339816, 242306, 340600,
                                                                       162704, 163244, 252302,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 354152, 0, 3,
                                                                       340600, 242894, 341384,
                                                                       163244, 163784, 253058,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 355160, 0, 3,
                                                                       341384, 243482, 342168,
                                                                       163784, 164324, 253814,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 356168, 0, 3,
                                                                       342168, 244070, 342952,
                                                                       164324, 164864, 254570,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 357176, 0, 3,
                                                                       342952, 244658, 343736,
                                                                       164864, 165404, 255326,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 358184, 0, 3,
                                                                       343736, 245246, 344520,
                                                                       165404, 165944, 256082,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 359192, 0, 3,
                                                                       344520, 245834, 345304,
                                                                       165944, 166484, 256838,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 360200, 0, 3,
                                                                       346088, 247010, 347096,
                                                                       167564, 168239, 257594,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 361460, 0, 3,
                                                                       347096, 247766, 348104,
                                                                       168239, 168914, 258539,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 362720, 0, 3,
                                                                       348104, 248522, 349112,
                                                                       168914, 169589, 259484,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 363980, 0, 3,
                                                                       349112, 249278, 350120,
                                                                       169589, 170264, 260429,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 365240, 0, 3,
                                                                       350120, 250034, 351128,
                                                                       170264, 170939, 261374,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 366500, 0, 3,
                                                                       351128, 250790, 352136,
                                                                       170939, 171614, 262319,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 367760, 0, 3,
                                                                       353144, 252302, 354152,
                                                                       172964, 173639, 263264,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 369020, 0, 3,
                                                                       354152, 253058, 355160,
                                                                       173639, 174314, 264209,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 370280, 0, 3,
                                                                       355160, 253814, 356168,
                                                                       174314, 174989, 265154,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 371540, 0, 3,
                                                                       356168, 254570, 357176,
                                                                       174989, 175664, 266099,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 372800, 0, 3,
                                                                       357176, 255326, 358184,
                                                                       175664, 176339, 267044,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 374060, 0, 3,
                                                                       358184, 256082, 359192,
                                                                       176339, 177014, 267989,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 375320, 0, 3,
                                                                       360200, 257594, 361460,
                                                                       178364, 179189, 268934,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 376860, 0, 3,
                                                                       361460, 258539, 362720,
                                                                       179189, 180014, 270089,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 378400, 0, 3,
                                                                       362720, 259484, 363980,
                                                                       180014, 180839, 271244,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 379940, 0, 3,
                                                                       363980, 260429, 365240,
                                                                       180839, 181664, 272399,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 381480, 0, 3,
                                                                       365240, 261374, 366500,
                                                                       181664, 182489, 273554,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 383020, 0, 3,
                                                                       367760, 263264, 369020,
                                                                       184139, 184964, 274709,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 384560, 0, 3,
                                                                       369020, 264209, 370280,
                                                                       184964, 185789, 275864,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 386100, 0, 3,
                                                                       370280, 265154, 371540,
                                                                       185789, 186614, 277019,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 387640, 0, 3,
                                                                       371540, 266099, 372800,
                                                                       186614, 187439, 278174,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 389180, 0, 3,
                                                                       372800, 267044, 374060,
                                                                       187439, 188264, 279329,
                                                                       ncols, gamma, p, q);

                    compute_prim_sni_three_center_electron_repulsion_0(buffer, 390720, 0, 3,
                                                                       375320, 268934, 376860,
                                                                       189914, 190904, 280484,
                                                                       ncols, gamma, p, q);

                    compute_prim_sni_three_center_electron_repulsion_0(buffer, 392568, 0, 3,
                                                                       376860, 270089, 378400,
                                                                       190904, 191894, 281870,
                                                                       ncols, gamma, p, q);

                    compute_prim_sni_three_center_electron_repulsion_0(buffer, 394416, 0, 3,
                                                                       378400, 271244, 379940,
                                                                       191894, 192884, 283256,
                                                                       ncols, gamma, p, q);

                    compute_prim_sni_three_center_electron_repulsion_0(buffer, 396264, 0, 3,
                                                                       379940, 272399, 381480,
                                                                       192884, 193874, 284642,
                                                                       ncols, gamma, p, q);

                    compute_prim_sni_three_center_electron_repulsion_0(buffer, 398112, 0, 3,
                                                                       383020, 274709, 384560,
                                                                       195854, 196844, 286028,
                                                                       ncols, gamma, p, q);

                    compute_prim_sni_three_center_electron_repulsion_0(buffer, 399960, 0, 3,
                                                                       384560, 275864, 386100,
                                                                       196844, 197834, 287414,
                                                                       ncols, gamma, p, q);

                    compute_prim_sni_three_center_electron_repulsion_0(buffer, 401808, 0, 3,
                                                                       386100, 277019, 387640,
                                                                       197834, 198824, 288800,
                                                                       ncols, gamma, p, q);

                    compute_prim_sni_three_center_electron_repulsion_0(buffer, 403656, 0, 3,
                                                                       387640, 278174, 389180,
                                                                       198824, 199814, 290186,
                                                                       ncols, gamma, p, q);

                    compute_prim_soi_three_center_electron_repulsion_0(buffer, 405504, 0, 3,
                                                                       390720, 280484, 392568,
                                                                       201794, 202964, 291572,
                                                                       ncols, gamma, p, q);

                    compute_prim_soi_three_center_electron_repulsion_0(buffer, 407688, 0, 3,
                                                                       392568, 281870, 394416,
                                                                       202964, 204134, 293210,
                                                                       ncols, gamma, p, q);

                    compute_prim_soi_three_center_electron_repulsion_0(buffer, 409872, 0, 3,
                                                                       394416, 283256, 396264,
                                                                       204134, 205304, 294848,
                                                                       ncols, gamma, p, q);

                    compute_prim_soi_three_center_electron_repulsion_0(buffer, 412056, 0, 3,
                                                                       398112, 286028, 399960,
                                                                       207644, 208814, 296486,
                                                                       ncols, gamma, p, q);

                    compute_prim_soi_three_center_electron_repulsion_0(buffer, 414240, 0, 3,
                                                                       399960, 287414, 401808,
                                                                       208814, 209984, 298124,
                                                                       ncols, gamma, p, q);

                    compute_prim_soi_three_center_electron_repulsion_0(buffer, 416424, 0, 3,
                                                                       401808, 288800, 403656,
                                                                       209984, 211154, 299762,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 418608, 3, 213494,
                                                                       213515, 301456, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 418644, 3, 213515,
                                                                       213536, 301484, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 418680, 3, 213536,
                                                                       213557, 301512, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 418716, 3, 213557,
                                                                       213578, 301540, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 418752, 3, 213578,
                                                                       213599, 301568, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 418788, 3, 213599,
                                                                       213620, 301596, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 418824, 3, 213620,
                                                                       213641, 301624, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 418860, 3, 213641,
                                                                       213662, 301652, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 418896, 3, 213662,
                                                                       213683, 301680, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 418932, 3, 213683,
                                                                       213704, 301708, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 418968, 3, 213704,
                                                                       213725, 301736, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 419004, 3, 213725,
                                                                       213746, 301764, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 419040, 3, 213788,
                                                                       213809, 301848, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 419076, 3, 213809,
                                                                       213830, 301876, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 419112, 3, 213830,
                                                                       213851, 301904, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 419148, 3, 213851,
                                                                       213872, 301932, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 419184, 3, 213872,
                                                                       213893, 301960, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 419220, 3, 213893,
                                                                       213914, 301988, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 419256, 3, 213914,
                                                                       213935, 302016, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 419292, 3, 213935,
                                                                       213956, 302044, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 419328, 3, 213956,
                                                                       213977, 302072, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 419364, 3, 213977,
                                                                       213998, 302100, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 419400, 3, 213998,
                                                                       214019, 302128, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 419436, 3, 214019,
                                                                       214040, 302156, ncols,
                                                                       gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 419472, 0, 3,
                                                                       418608, 301456, 418644,
                                                                       214082, 214145, 302352,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 419580, 0, 3,
                                                                       418644, 301484, 418680,
                                                                       214145, 214208, 302436,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 419688, 0, 3,
                                                                       418680, 301512, 418716,
                                                                       214208, 214271, 302520,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 419796, 0, 3,
                                                                       418716, 301540, 418752,
                                                                       214271, 214334, 302604,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 419904, 0, 3,
                                                                       418752, 301568, 418788,
                                                                       214334, 214397, 302688,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 420012, 0, 3,
                                                                       418788, 301596, 418824,
                                                                       214397, 214460, 302772,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 420120, 0, 3,
                                                                       418824, 301624, 418860,
                                                                       214460, 214523, 302856,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 420228, 0, 3,
                                                                       418860, 301652, 418896,
                                                                       214523, 214586, 302940,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 420336, 0, 3,
                                                                       418896, 301680, 418932,
                                                                       214586, 214649, 303024,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 420444, 0, 3,
                                                                       418932, 301708, 418968,
                                                                       214649, 214712, 303108,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 420552, 0, 3,
                                                                       418968, 301736, 419004,
                                                                       214712, 214775, 303192,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 420660, 0, 3,
                                                                       419040, 301848, 419076,
                                                                       214901, 214964, 303444,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 420768, 0, 3,
                                                                       419076, 301876, 419112,
                                                                       214964, 215027, 303528,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 420876, 0, 3,
                                                                       419112, 301904, 419148,
                                                                       215027, 215090, 303612,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 420984, 0, 3,
                                                                       419148, 301932, 419184,
                                                                       215090, 215153, 303696,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 421092, 0, 3,
                                                                       419184, 301960, 419220,
                                                                       215153, 215216, 303780,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 421200, 0, 3,
                                                                       419220, 301988, 419256,
                                                                       215216, 215279, 303864,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 421308, 0, 3,
                                                                       419256, 302016, 419292,
                                                                       215279, 215342, 303948,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 421416, 0, 3,
                                                                       419292, 302044, 419328,
                                                                       215342, 215405, 304032,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 421524, 0, 3,
                                                                       419328, 302072, 419364,
                                                                       215405, 215468, 304116,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 421632, 0, 3,
                                                                       419364, 302100, 419400,
                                                                       215468, 215531, 304200,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 421740, 0, 3,
                                                                       419400, 302128, 419436,
                                                                       215531, 215594, 304284,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 421848, 0, 3,
                                                                       419472, 302352, 419580,
                                                                       215720, 215846, 304704,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 422064, 0, 3,
                                                                       419580, 302436, 419688,
                                                                       215846, 215972, 304872,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 422280, 0, 3,
                                                                       419688, 302520, 419796,
                                                                       215972, 216098, 305040,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 422496, 0, 3,
                                                                       419796, 302604, 419904,
                                                                       216098, 216224, 305208,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 422712, 0, 3,
                                                                       419904, 302688, 420012,
                                                                       216224, 216350, 305376,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 422928, 0, 3,
                                                                       420012, 302772, 420120,
                                                                       216350, 216476, 305544,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 423144, 0, 3,
                                                                       420120, 302856, 420228,
                                                                       216476, 216602, 305712,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 423360, 0, 3,
                                                                       420228, 302940, 420336,
                                                                       216602, 216728, 305880,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 423576, 0, 3,
                                                                       420336, 303024, 420444,
                                                                       216728, 216854, 306048,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 423792, 0, 3,
                                                                       420444, 303108, 420552,
                                                                       216854, 216980, 306216,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 424008, 0, 3,
                                                                       420660, 303444, 420768,
                                                                       217232, 217358, 306720,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 424224, 0, 3,
                                                                       420768, 303528, 420876,
                                                                       217358, 217484, 306888,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 424440, 0, 3,
                                                                       420876, 303612, 420984,
                                                                       217484, 217610, 307056,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 424656, 0, 3,
                                                                       420984, 303696, 421092,
                                                                       217610, 217736, 307224,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 424872, 0, 3,
                                                                       421092, 303780, 421200,
                                                                       217736, 217862, 307392,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 425088, 0, 3,
                                                                       421200, 303864, 421308,
                                                                       217862, 217988, 307560,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 425304, 0, 3,
                                                                       421308, 303948, 421416,
                                                                       217988, 218114, 307728,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 425520, 0, 3,
                                                                       421416, 304032, 421524,
                                                                       218114, 218240, 307896,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 425736, 0, 3,
                                                                       421524, 304116, 421632,
                                                                       218240, 218366, 308064,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 425952, 0, 3,
                                                                       421632, 304200, 421740,
                                                                       218366, 218492, 308232,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 426168, 0, 3,
                                                                       421848, 304704, 422064,
                                                                       218744, 218954, 308960,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 426528, 0, 3,
                                                                       422064, 304872, 422280,
                                                                       218954, 219164, 309240,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 426888, 0, 3,
                                                                       422280, 305040, 422496,
                                                                       219164, 219374, 309520,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 427248, 0, 3,
                                                                       422496, 305208, 422712,
                                                                       219374, 219584, 309800,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 427608, 0, 3,
                                                                       422712, 305376, 422928,
                                                                       219584, 219794, 310080,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 427968, 0, 3,
                                                                       422928, 305544, 423144,
                                                                       219794, 220004, 310360,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 428328, 0, 3,
                                                                       423144, 305712, 423360,
                                                                       220004, 220214, 310640,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 428688, 0, 3,
                                                                       423360, 305880, 423576,
                                                                       220214, 220424, 310920,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 429048, 0, 3,
                                                                       423576, 306048, 423792,
                                                                       220424, 220634, 311200,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 429408, 0, 3,
                                                                       424008, 306720, 424224,
                                                                       221054, 221264, 312040,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 429768, 0, 3,
                                                                       424224, 306888, 424440,
                                                                       221264, 221474, 312320,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 430128, 0, 3,
                                                                       424440, 307056, 424656,
                                                                       221474, 221684, 312600,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 430488, 0, 3,
                                                                       424656, 307224, 424872,
                                                                       221684, 221894, 312880,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 430848, 0, 3,
                                                                       424872, 307392, 425088,
                                                                       221894, 222104, 313160,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 431208, 0, 3,
                                                                       425088, 307560, 425304,
                                                                       222104, 222314, 313440,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 431568, 0, 3,
                                                                       425304, 307728, 425520,
                                                                       222314, 222524, 313720,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 431928, 0, 3,
                                                                       425520, 307896, 425736,
                                                                       222524, 222734, 314000,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 432288, 0, 3,
                                                                       425736, 308064, 425952,
                                                                       222734, 222944, 314280,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 432648, 0, 3,
                                                                       426168, 308960, 426528,
                                                                       223364, 223679, 315400,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 433188, 0, 3,
                                                                       426528, 309240, 426888,
                                                                       223679, 223994, 315820,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 433728, 0, 3,
                                                                       426888, 309520, 427248,
                                                                       223994, 224309, 316240,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 434268, 0, 3,
                                                                       427248, 309800, 427608,
                                                                       224309, 224624, 316660,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 434808, 0, 3,
                                                                       427608, 310080, 427968,
                                                                       224624, 224939, 317080,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 435348, 0, 3,
                                                                       427968, 310360, 428328,
                                                                       224939, 225254, 317500,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 435888, 0, 3,
                                                                       428328, 310640, 428688,
                                                                       225254, 225569, 317920,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 436428, 0, 3,
                                                                       428688, 310920, 429048,
                                                                       225569, 225884, 318340,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 436968, 0, 3,
                                                                       429408, 312040, 429768,
                                                                       226514, 226829, 319600,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 437508, 0, 3,
                                                                       429768, 312320, 430128,
                                                                       226829, 227144, 320020,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 438048, 0, 3,
                                                                       430128, 312600, 430488,
                                                                       227144, 227459, 320440,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 438588, 0, 3,
                                                                       430488, 312880, 430848,
                                                                       227459, 227774, 320860,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 439128, 0, 3,
                                                                       430848, 313160, 431208,
                                                                       227774, 228089, 321280,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 439668, 0, 3,
                                                                       431208, 313440, 431568,
                                                                       228089, 228404, 321700,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 440208, 0, 3,
                                                                       431568, 313720, 431928,
                                                                       228404, 228719, 322120,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 440748, 0, 3,
                                                                       431928, 314000, 432288,
                                                                       228719, 229034, 322540,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 441288, 0, 3,
                                                                       432648, 315400, 433188,
                                                                       229664, 230105, 324136,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 442044, 0, 3,
                                                                       433188, 315820, 433728,
                                                                       230105, 230546, 324724,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 442800, 0, 3,
                                                                       433728, 316240, 434268,
                                                                       230546, 230987, 325312,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 443556, 0, 3,
                                                                       434268, 316660, 434808,
                                                                       230987, 231428, 325900,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 444312, 0, 3,
                                                                       434808, 317080, 435348,
                                                                       231428, 231869, 326488,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 445068, 0, 3,
                                                                       435348, 317500, 435888,
                                                                       231869, 232310, 327076,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 445824, 0, 3,
                                                                       435888, 317920, 436428,
                                                                       232310, 232751, 327664,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 446580, 0, 3,
                                                                       436968, 319600, 437508,
                                                                       233633, 234074, 329428,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 447336, 0, 3,
                                                                       437508, 320020, 438048,
                                                                       234074, 234515, 330016,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 448092, 0, 3,
                                                                       438048, 320440, 438588,
                                                                       234515, 234956, 330604,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 448848, 0, 3,
                                                                       438588, 320860, 439128,
                                                                       234956, 235397, 331192,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 449604, 0, 3,
                                                                       439128, 321280, 439668,
                                                                       235397, 235838, 331780,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 450360, 0, 3,
                                                                       439668, 321700, 440208,
                                                                       235838, 236279, 332368,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 451116, 0, 3,
                                                                       440208, 322120, 440748,
                                                                       236279, 236720, 332956,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 451872, 0, 3,
                                                                       441288, 324136, 442044,
                                                                       237602, 238190, 335112,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 452880, 0, 3,
                                                                       442044, 324724, 442800,
                                                                       238190, 238778, 335896,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 453888, 0, 3,
                                                                       442800, 325312, 443556,
                                                                       238778, 239366, 336680,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 454896, 0, 3,
                                                                       443556, 325900, 444312,
                                                                       239366, 239954, 337464,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 455904, 0, 3,
                                                                       444312, 326488, 445068,
                                                                       239954, 240542, 338248,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 456912, 0, 3,
                                                                       445068, 327076, 445824,
                                                                       240542, 241130, 339032,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 457920, 0, 3,
                                                                       446580, 329428, 447336,
                                                                       242306, 242894, 341384,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 458928, 0, 3,
                                                                       447336, 330016, 448092,
                                                                       242894, 243482, 342168,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 459936, 0, 3,
                                                                       448092, 330604, 448848,
                                                                       243482, 244070, 342952,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 460944, 0, 3,
                                                                       448848, 331192, 449604,
                                                                       244070, 244658, 343736,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 461952, 0, 3,
                                                                       449604, 331780, 450360,
                                                                       244658, 245246, 344520,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 462960, 0, 3,
                                                                       450360, 332368, 451116,
                                                                       245246, 245834, 345304,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 463968, 0, 3,
                                                                       451872, 335112, 452880,
                                                                       247010, 247766, 348104,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 465264, 0, 3,
                                                                       452880, 335896, 453888,
                                                                       247766, 248522, 349112,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 466560, 0, 3,
                                                                       453888, 336680, 454896,
                                                                       248522, 249278, 350120,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 467856, 0, 3,
                                                                       454896, 337464, 455904,
                                                                       249278, 250034, 351128,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 469152, 0, 3,
                                                                       455904, 338248, 456912,
                                                                       250034, 250790, 352136,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 470448, 0, 3,
                                                                       457920, 341384, 458928,
                                                                       252302, 253058, 355160,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 471744, 0, 3,
                                                                       458928, 342168, 459936,
                                                                       253058, 253814, 356168,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 473040, 0, 3,
                                                                       459936, 342952, 460944,
                                                                       253814, 254570, 357176,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 474336, 0, 3,
                                                                       460944, 343736, 461952,
                                                                       254570, 255326, 358184,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 475632, 0, 3,
                                                                       461952, 344520, 462960,
                                                                       255326, 256082, 359192,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 476928, 0, 3,
                                                                       463968, 348104, 465264,
                                                                       257594, 258539, 362720,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 478548, 0, 3,
                                                                       465264, 349112, 466560,
                                                                       258539, 259484, 363980,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 480168, 0, 3,
                                                                       466560, 350120, 467856,
                                                                       259484, 260429, 365240,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 481788, 0, 3,
                                                                       467856, 351128, 469152,
                                                                       260429, 261374, 366500,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 483408, 0, 3,
                                                                       470448, 355160, 471744,
                                                                       263264, 264209, 370280,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 485028, 0, 3,
                                                                       471744, 356168, 473040,
                                                                       264209, 265154, 371540,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 486648, 0, 3,
                                                                       473040, 357176, 474336,
                                                                       265154, 266099, 372800,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 488268, 0, 3,
                                                                       474336, 358184, 475632,
                                                                       266099, 267044, 374060,
                                                                       ncols, gamma, p, q);

                    compute_prim_smk_three_center_electron_repulsion_0(buffer, 489888, 0, 3,
                                                                       476928, 362720, 478548,
                                                                       268934, 270089, 378400,
                                                                       ncols, gamma, p, q);

                    compute_prim_smk_three_center_electron_repulsion_0(buffer, 491868, 0, 3,
                                                                       478548, 363980, 480168,
                                                                       270089, 271244, 379940,
                                                                       ncols, gamma, p, q);

                    compute_prim_smk_three_center_electron_repulsion_0(buffer, 493848, 0, 3,
                                                                       480168, 365240, 481788,
                                                                       271244, 272399, 381480,
                                                                       ncols, gamma, p, q);

                    compute_prim_smk_three_center_electron_repulsion_0(buffer, 495828, 0, 3,
                                                                       483408, 370280, 485028,
                                                                       274709, 275864, 386100,
                                                                       ncols, gamma, p, q);

                    compute_prim_smk_three_center_electron_repulsion_0(buffer, 497808, 0, 3,
                                                                       485028, 371540, 486648,
                                                                       275864, 277019, 387640,
                                                                       ncols, gamma, p, q);

                    compute_prim_smk_three_center_electron_repulsion_0(buffer, 499788, 0, 3,
                                                                       486648, 372800, 488268,
                                                                       277019, 278174, 389180,
                                                                       ncols, gamma, p, q);

                    compute_prim_snk_three_center_electron_repulsion_0(buffer, 501768, 0, 3,
                                                                       489888, 378400, 491868,
                                                                       280484, 281870, 394416,
                                                                       ncols, gamma, p, q);

                    compute_prim_snk_three_center_electron_repulsion_0(buffer, 504144, 0, 3,
                                                                       491868, 379940, 493848,
                                                                       281870, 283256, 396264,
                                                                       ncols, gamma, p, q);

                    compute_prim_snk_three_center_electron_repulsion_0(buffer, 506520, 0, 3,
                                                                       495828, 386100, 497808,
                                                                       286028, 287414, 401808,
                                                                       ncols, gamma, p, q);

                    compute_prim_snk_three_center_electron_repulsion_0(buffer, 508896, 0, 3,
                                                                       497808, 387640, 499788,
                                                                       287414, 288800, 403656,
                                                                       ncols, gamma, p, q);

                    compute_prim_sok_three_center_electron_repulsion_0(buffer, 511272, 0, 3,
                                                                       501768, 394416, 504144,
                                                                       291572, 293210, 409872,
                                                                       ncols, gamma, p, q);

                    compute_prim_sok_three_center_electron_repulsion_0(buffer, 514080, 0, 3,
                                                                       506520, 401808, 508896,
                                                                       296486, 298124, 416424,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 516888, 3, 301400,
                                                                       301428, 418608, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 516933, 3, 301428,
                                                                       301456, 418644, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 516978, 3, 301456,
                                                                       301484, 418680, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 517023, 3, 301484,
                                                                       301512, 418716, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 517068, 3, 301512,
                                                                       301540, 418752, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 517113, 3, 301540,
                                                                       301568, 418788, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 517158, 3, 301568,
                                                                       301596, 418824, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 517203, 3, 301596,
                                                                       301624, 418860, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 517248, 3, 301624,
                                                                       301652, 418896, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 517293, 3, 301652,
                                                                       301680, 418932, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 517338, 3, 301680,
                                                                       301708, 418968, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 517383, 3, 301708,
                                                                       301736, 419004, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 517428, 3, 301792,
                                                                       301820, 419040, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 517473, 3, 301820,
                                                                       301848, 419076, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 517518, 3, 301848,
                                                                       301876, 419112, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 517563, 3, 301876,
                                                                       301904, 419148, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 517608, 3, 301904,
                                                                       301932, 419184, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 517653, 3, 301932,
                                                                       301960, 419220, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 517698, 3, 301960,
                                                                       301988, 419256, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 517743, 3, 301988,
                                                                       302016, 419292, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 517788, 3, 302016,
                                                                       302044, 419328, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 517833, 3, 302044,
                                                                       302072, 419364, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 517878, 3, 302072,
                                                                       302100, 419400, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 517923, 3, 302100,
                                                                       302128, 419436, ncols,
                                                                       gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 517968, 0, 3,
                                                                       516888, 418608, 516933,
                                                                       302184, 302268, 419472,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 518103, 0, 3,
                                                                       516933, 418644, 516978,
                                                                       302268, 302352, 419580,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 518238, 0, 3,
                                                                       516978, 418680, 517023,
                                                                       302352, 302436, 419688,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 518373, 0, 3,
                                                                       517023, 418716, 517068,
                                                                       302436, 302520, 419796,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 518508, 0, 3,
                                                                       517068, 418752, 517113,
                                                                       302520, 302604, 419904,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 518643, 0, 3,
                                                                       517113, 418788, 517158,
                                                                       302604, 302688, 420012,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 518778, 0, 3,
                                                                       517158, 418824, 517203,
                                                                       302688, 302772, 420120,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 518913, 0, 3,
                                                                       517203, 418860, 517248,
                                                                       302772, 302856, 420228,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 519048, 0, 3,
                                                                       517248, 418896, 517293,
                                                                       302856, 302940, 420336,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 519183, 0, 3,
                                                                       517293, 418932, 517338,
                                                                       302940, 303024, 420444,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 519318, 0, 3,
                                                                       517338, 418968, 517383,
                                                                       303024, 303108, 420552,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 519453, 0, 3,
                                                                       517428, 419040, 517473,
                                                                       303276, 303360, 420660,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 519588, 0, 3,
                                                                       517473, 419076, 517518,
                                                                       303360, 303444, 420768,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 519723, 0, 3,
                                                                       517518, 419112, 517563,
                                                                       303444, 303528, 420876,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 519858, 0, 3,
                                                                       517563, 419148, 517608,
                                                                       303528, 303612, 420984,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 519993, 0, 3,
                                                                       517608, 419184, 517653,
                                                                       303612, 303696, 421092,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 520128, 0, 3,
                                                                       517653, 419220, 517698,
                                                                       303696, 303780, 421200,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 520263, 0, 3,
                                                                       517698, 419256, 517743,
                                                                       303780, 303864, 421308,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 520398, 0, 3,
                                                                       517743, 419292, 517788,
                                                                       303864, 303948, 421416,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 520533, 0, 3,
                                                                       517788, 419328, 517833,
                                                                       303948, 304032, 421524,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 520668, 0, 3,
                                                                       517833, 419364, 517878,
                                                                       304032, 304116, 421632,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 520803, 0, 3,
                                                                       517878, 419400, 517923,
                                                                       304116, 304200, 421740,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 520938, 0, 3,
                                                                       517968, 419472, 518103,
                                                                       304368, 304536, 421848,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 521208, 0, 3,
                                                                       518103, 419580, 518238,
                                                                       304536, 304704, 422064,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 521478, 0, 3,
                                                                       518238, 419688, 518373,
                                                                       304704, 304872, 422280,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 521748, 0, 3,
                                                                       518373, 419796, 518508,
                                                                       304872, 305040, 422496,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 522018, 0, 3,
                                                                       518508, 419904, 518643,
                                                                       305040, 305208, 422712,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 522288, 0, 3,
                                                                       518643, 420012, 518778,
                                                                       305208, 305376, 422928,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 522558, 0, 3,
                                                                       518778, 420120, 518913,
                                                                       305376, 305544, 423144,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 522828, 0, 3,
                                                                       518913, 420228, 519048,
                                                                       305544, 305712, 423360,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 523098, 0, 3,
                                                                       519048, 420336, 519183,
                                                                       305712, 305880, 423576,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 523368, 0, 3,
                                                                       519183, 420444, 519318,
                                                                       305880, 306048, 423792,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 523638, 0, 3,
                                                                       519453, 420660, 519588,
                                                                       306384, 306552, 424008,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 523908, 0, 3,
                                                                       519588, 420768, 519723,
                                                                       306552, 306720, 424224,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 524178, 0, 3,
                                                                       519723, 420876, 519858,
                                                                       306720, 306888, 424440,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 524448, 0, 3,
                                                                       519858, 420984, 519993,
                                                                       306888, 307056, 424656,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 524718, 0, 3,
                                                                       519993, 421092, 520128,
                                                                       307056, 307224, 424872,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 524988, 0, 3,
                                                                       520128, 421200, 520263,
                                                                       307224, 307392, 425088,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 525258, 0, 3,
                                                                       520263, 421308, 520398,
                                                                       307392, 307560, 425304,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 525528, 0, 3,
                                                                       520398, 421416, 520533,
                                                                       307560, 307728, 425520,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 525798, 0, 3,
                                                                       520533, 421524, 520668,
                                                                       307728, 307896, 425736,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 526068, 0, 3,
                                                                       520668, 421632, 520803,
                                                                       307896, 308064, 425952,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 526338, 0, 3,
                                                                       520938, 421848, 521208,
                                                                       308400, 308680, 426168,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 526788, 0, 3,
                                                                       521208, 422064, 521478,
                                                                       308680, 308960, 426528,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 527238, 0, 3,
                                                                       521478, 422280, 521748,
                                                                       308960, 309240, 426888,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 527688, 0, 3,
                                                                       521748, 422496, 522018,
                                                                       309240, 309520, 427248,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 528138, 0, 3,
                                                                       522018, 422712, 522288,
                                                                       309520, 309800, 427608,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 528588, 0, 3,
                                                                       522288, 422928, 522558,
                                                                       309800, 310080, 427968,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 529038, 0, 3,
                                                                       522558, 423144, 522828,
                                                                       310080, 310360, 428328,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 529488, 0, 3,
                                                                       522828, 423360, 523098,
                                                                       310360, 310640, 428688,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 529938, 0, 3,
                                                                       523098, 423576, 523368,
                                                                       310640, 310920, 429048,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 530388, 0, 3,
                                                                       523638, 424008, 523908,
                                                                       311480, 311760, 429408,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 530838, 0, 3,
                                                                       523908, 424224, 524178,
                                                                       311760, 312040, 429768,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 531288, 0, 3,
                                                                       524178, 424440, 524448,
                                                                       312040, 312320, 430128,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 531738, 0, 3,
                                                                       524448, 424656, 524718,
                                                                       312320, 312600, 430488,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 532188, 0, 3,
                                                                       524718, 424872, 524988,
                                                                       312600, 312880, 430848,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 532638, 0, 3,
                                                                       524988, 425088, 525258,
                                                                       312880, 313160, 431208,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 533088, 0, 3,
                                                                       525258, 425304, 525528,
                                                                       313160, 313440, 431568,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 533538, 0, 3,
                                                                       525528, 425520, 525798,
                                                                       313440, 313720, 431928,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 533988, 0, 3,
                                                                       525798, 425736, 526068,
                                                                       313720, 314000, 432288,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 534438, 0, 3,
                                                                       526338, 426168, 526788,
                                                                       314560, 314980, 432648,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 535113, 0, 3,
                                                                       526788, 426528, 527238,
                                                                       314980, 315400, 433188,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 535788, 0, 3,
                                                                       527238, 426888, 527688,
                                                                       315400, 315820, 433728,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 536463, 0, 3,
                                                                       527688, 427248, 528138,
                                                                       315820, 316240, 434268,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 537138, 0, 3,
                                                                       528138, 427608, 528588,
                                                                       316240, 316660, 434808,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 537813, 0, 3,
                                                                       528588, 427968, 529038,
                                                                       316660, 317080, 435348,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 538488, 0, 3,
                                                                       529038, 428328, 529488,
                                                                       317080, 317500, 435888,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 539163, 0, 3,
                                                                       529488, 428688, 529938,
                                                                       317500, 317920, 436428,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 539838, 0, 3,
                                                                       530388, 429408, 530838,
                                                                       318760, 319180, 436968,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 540513, 0, 3,
                                                                       530838, 429768, 531288,
                                                                       319180, 319600, 437508,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 541188, 0, 3,
                                                                       531288, 430128, 531738,
                                                                       319600, 320020, 438048,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 541863, 0, 3,
                                                                       531738, 430488, 532188,
                                                                       320020, 320440, 438588,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 542538, 0, 3,
                                                                       532188, 430848, 532638,
                                                                       320440, 320860, 439128,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 543213, 0, 3,
                                                                       532638, 431208, 533088,
                                                                       320860, 321280, 439668,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 543888, 0, 3,
                                                                       533088, 431568, 533538,
                                                                       321280, 321700, 440208,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 544563, 0, 3,
                                                                       533538, 431928, 533988,
                                                                       321700, 322120, 440748,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 545238, 0, 3,
                                                                       534438, 432648, 535113,
                                                                       322960, 323548, 441288,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 546183, 0, 3,
                                                                       535113, 433188, 535788,
                                                                       323548, 324136, 442044,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 547128, 0, 3,
                                                                       535788, 433728, 536463,
                                                                       324136, 324724, 442800,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 548073, 0, 3,
                                                                       536463, 434268, 537138,
                                                                       324724, 325312, 443556,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 549018, 0, 3,
                                                                       537138, 434808, 537813,
                                                                       325312, 325900, 444312,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 549963, 0, 3,
                                                                       537813, 435348, 538488,
                                                                       325900, 326488, 445068,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 550908, 0, 3,
                                                                       538488, 435888, 539163,
                                                                       326488, 327076, 445824,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 551853, 0, 3,
                                                                       539838, 436968, 540513,
                                                                       328252, 328840, 446580,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 552798, 0, 3,
                                                                       540513, 437508, 541188,
                                                                       328840, 329428, 447336,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 553743, 0, 3,
                                                                       541188, 438048, 541863,
                                                                       329428, 330016, 448092,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 554688, 0, 3,
                                                                       541863, 438588, 542538,
                                                                       330016, 330604, 448848,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 555633, 0, 3,
                                                                       542538, 439128, 543213,
                                                                       330604, 331192, 449604,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 556578, 0, 3,
                                                                       543213, 439668, 543888,
                                                                       331192, 331780, 450360,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 557523, 0, 3,
                                                                       543888, 440208, 544563,
                                                                       331780, 332368, 451116,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 558468, 0, 3,
                                                                       545238, 441288, 546183,
                                                                       333544, 334328, 451872,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 559728, 0, 3,
                                                                       546183, 442044, 547128,
                                                                       334328, 335112, 452880,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 560988, 0, 3,
                                                                       547128, 442800, 548073,
                                                                       335112, 335896, 453888,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 562248, 0, 3,
                                                                       548073, 443556, 549018,
                                                                       335896, 336680, 454896,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 563508, 0, 3,
                                                                       549018, 444312, 549963,
                                                                       336680, 337464, 455904,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 564768, 0, 3,
                                                                       549963, 445068, 550908,
                                                                       337464, 338248, 456912,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 566028, 0, 3,
                                                                       551853, 446580, 552798,
                                                                       339816, 340600, 457920,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 567288, 0, 3,
                                                                       552798, 447336, 553743,
                                                                       340600, 341384, 458928,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 568548, 0, 3,
                                                                       553743, 448092, 554688,
                                                                       341384, 342168, 459936,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 569808, 0, 3,
                                                                       554688, 448848, 555633,
                                                                       342168, 342952, 460944,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 571068, 0, 3,
                                                                       555633, 449604, 556578,
                                                                       342952, 343736, 461952,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 572328, 0, 3,
                                                                       556578, 450360, 557523,
                                                                       343736, 344520, 462960,
                                                                       ncols, gamma, p, q);

                    compute_prim_skl_three_center_electron_repulsion_0(buffer, 573588, 0, 3,
                                                                       558468, 451872, 559728,
                                                                       346088, 347096, 463968,
                                                                       ncols, gamma, p, q);

                    compute_prim_skl_three_center_electron_repulsion_0(buffer, 575208, 0, 3,
                                                                       559728, 452880, 560988,
                                                                       347096, 348104, 465264,
                                                                       ncols, gamma, p, q);

                    compute_prim_skl_three_center_electron_repulsion_0(buffer, 576828, 0, 3,
                                                                       560988, 453888, 562248,
                                                                       348104, 349112, 466560,
                                                                       ncols, gamma, p, q);

                    compute_prim_skl_three_center_electron_repulsion_0(buffer, 578448, 0, 3,
                                                                       562248, 454896, 563508,
                                                                       349112, 350120, 467856,
                                                                       ncols, gamma, p, q);

                    compute_prim_skl_three_center_electron_repulsion_0(buffer, 580068, 0, 3,
                                                                       563508, 455904, 564768,
                                                                       350120, 351128, 469152,
                                                                       ncols, gamma, p, q);

                    compute_prim_skl_three_center_electron_repulsion_0(buffer, 581688, 0, 3,
                                                                       566028, 457920, 567288,
                                                                       353144, 354152, 470448,
                                                                       ncols, gamma, p, q);

                    compute_prim_skl_three_center_electron_repulsion_0(buffer, 583308, 0, 3,
                                                                       567288, 458928, 568548,
                                                                       354152, 355160, 471744,
                                                                       ncols, gamma, p, q);

                    compute_prim_skl_three_center_electron_repulsion_0(buffer, 584928, 0, 3,
                                                                       568548, 459936, 569808,
                                                                       355160, 356168, 473040,
                                                                       ncols, gamma, p, q);

                    compute_prim_skl_three_center_electron_repulsion_0(buffer, 586548, 0, 3,
                                                                       569808, 460944, 571068,
                                                                       356168, 357176, 474336,
                                                                       ncols, gamma, p, q);

                    compute_prim_skl_three_center_electron_repulsion_0(buffer, 588168, 0, 3,
                                                                       571068, 461952, 572328,
                                                                       357176, 358184, 475632,
                                                                       ncols, gamma, p, q);

                    compute_prim_sll_three_center_electron_repulsion_0(buffer, 589788, 0, 3,
                                                                       573588, 463968, 575208,
                                                                       360200, 361460, 476928,
                                                                       ncols, gamma, p, q);

                    compute_prim_sll_three_center_electron_repulsion_0(buffer, 591813, 0, 3,
                                                                       575208, 465264, 576828,
                                                                       361460, 362720, 478548,
                                                                       ncols, gamma, p, q);

                    compute_prim_sll_three_center_electron_repulsion_0(buffer, 593838, 0, 3,
                                                                       576828, 466560, 578448,
                                                                       362720, 363980, 480168,
                                                                       ncols, gamma, p, q);

                    compute_prim_sll_three_center_electron_repulsion_0(buffer, 595863, 0, 3,
                                                                       578448, 467856, 580068,
                                                                       363980, 365240, 481788,
                                                                       ncols, gamma, p, q);

                    compute_prim_sll_three_center_electron_repulsion_0(buffer, 597888, 0, 3,
                                                                       581688, 470448, 583308,
                                                                       367760, 369020, 483408,
                                                                       ncols, gamma, p, q);

                    compute_prim_sll_three_center_electron_repulsion_0(buffer, 599913, 0, 3,
                                                                       583308, 471744, 584928,
                                                                       369020, 370280, 485028,
                                                                       ncols, gamma, p, q);

                    compute_prim_sll_three_center_electron_repulsion_0(buffer, 601938, 0, 3,
                                                                       584928, 473040, 586548,
                                                                       370280, 371540, 486648,
                                                                       ncols, gamma, p, q);

                    compute_prim_sll_three_center_electron_repulsion_0(buffer, 603963, 0, 3,
                                                                       586548, 474336, 588168,
                                                                       371540, 372800, 488268,
                                                                       ncols, gamma, p, q);

                    compute_prim_sml_three_center_electron_repulsion_0(buffer, 605988, 0, 3,
                                                                       589788, 476928, 591813,
                                                                       375320, 376860, 489888,
                                                                       ncols, gamma, p, q);

                    compute_prim_sml_three_center_electron_repulsion_0(buffer, 608463, 0, 3,
                                                                       591813, 478548, 593838,
                                                                       376860, 378400, 491868,
                                                                       ncols, gamma, p, q);

                    compute_prim_sml_three_center_electron_repulsion_0(buffer, 610938, 0, 3,
                                                                       593838, 480168, 595863,
                                                                       378400, 379940, 493848,
                                                                       ncols, gamma, p, q);

                    compute_prim_sml_three_center_electron_repulsion_0(buffer, 613413, 0, 3,
                                                                       597888, 483408, 599913,
                                                                       383020, 384560, 495828,
                                                                       ncols, gamma, p, q);

                    compute_prim_sml_three_center_electron_repulsion_0(buffer, 615888, 0, 3,
                                                                       599913, 485028, 601938,
                                                                       384560, 386100, 497808,
                                                                       ncols, gamma, p, q);

                    compute_prim_sml_three_center_electron_repulsion_0(buffer, 618363, 0, 3,
                                                                       601938, 486648, 603963,
                                                                       386100, 387640, 499788,
                                                                       ncols, gamma, p, q);

                    compute_prim_snl_three_center_electron_repulsion_0(buffer, 620838, 0, 3,
                                                                       605988, 489888, 608463,
                                                                       390720, 392568, 501768,
                                                                       ncols, gamma, p, q);

                    compute_prim_snl_three_center_electron_repulsion_0(buffer, 623808, 0, 3,
                                                                       608463, 491868, 610938,
                                                                       392568, 394416, 504144,
                                                                       ncols, gamma, p, q);

                    compute_prim_snl_three_center_electron_repulsion_0(buffer, 626778, 0, 3,
                                                                       613413, 495828, 615888,
                                                                       398112, 399960, 506520,
                                                                       ncols, gamma, p, q);

                    compute_prim_snl_three_center_electron_repulsion_0(buffer, 629748, 0, 3,
                                                                       615888, 497808, 618363,
                                                                       399960, 401808, 508896,
                                                                       ncols, gamma, p, q);

                    compute_prim_sol_three_center_electron_repulsion_0(buffer, 632718, 0, 3,
                                                                       620838, 501768, 623808,
                                                                       405504, 407688, 511272,
                                                                       ncols, gamma, p, q);

                    compute_prim_sol_three_center_electron_repulsion_0(buffer, 636228, 0, 3,
                                                                       626778, 506520, 629748,
                                                                       412056, 414240, 514080,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 639738, 558468, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 641474, 566028, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 643210, 573588, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 645442, 581688, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 647674, 589788, 2025, ncols);

                    simdfunc::contract_primitives(buffer, 650464, 597888, 2025, ncols);

                    simdfunc::contract_primitives(buffer, 653254, 605988, 2475, ncols);

                    simdfunc::contract_primitives(buffer, 656664, 613413, 2475, ncols);

                    simdfunc::contract_primitives(buffer, 660074, 620838, 2970, ncols);

                    simdfunc::contract_primitives(buffer, 664166, 626778, 2970, ncols);

                    simdfunc::contract_primitives(buffer, 668258, 632718, 3510, ncols);

                    simdfunc::contract_primitives(buffer, 673094, 636228, 3510, ncols);
                }
            }
        }

        simdtrf::transform_l_inner(buffer, 640998, 639738, 28, 1, nmax);

        simdtrf::transform_l_inner(buffer, 642734, 641474, 28, 1, nmax);

        simdtrf::transform_l_inner(buffer, 644830, 643210, 36, 1, nmax);

        simdtrf::transform_l_inner(buffer, 647062, 645442, 36, 1, nmax);

        simdtrf::transform_l_inner(buffer, 649699, 647674, 45, 1, nmax);

        simdtrf::transform_l_inner(buffer, 652489, 650464, 45, 1, nmax);

        simdtrf::transform_l_inner(buffer, 655729, 653254, 55, 1, nmax);

        simdtrf::transform_l_inner(buffer, 659139, 656664, 55, 1, nmax);

        simdtrf::transform_l_inner(buffer, 663044, 660074, 66, 1, nmax);

        simdtrf::transform_l_inner(buffer, 667136, 664166, 66, 1, nmax);

        simdtrf::transform_l_inner(buffer, 671768, 668258, 78, 1, nmax);

        simdtrf::transform_l_inner(buffer, 676604, 673094, 78, 1, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 677930, 640998, 644830, 17, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 679358, 642734, 647062, 17, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 680786, 644830, 649699, 17, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 682622, 647062, 652489, 17, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 684458, 649699, 655729, 17, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 686753, 652489, 659139, 17, nmax);

        simdtrf::compute_hrr_pm(buffer, coordinates, 689048, 655729, 663044, 17, nmax);

        simdtrf::compute_hrr_pm(buffer, coordinates, 691853, 659139, 667136, 17, nmax);

        simdtrf::compute_hrr_pn(buffer, coordinates, 694658, 663044, 671768, 17, nmax);

        simdtrf::compute_hrr_pn(buffer, coordinates, 698024, 667136, 676604, 17, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 701390, 677930, 680786, 17, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 704246, 679358, 682622, 17, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 707102, 680786, 684458, 17, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 710774, 682622, 686753, 17, nmax);

        simdtrf::compute_hrr_dl(buffer, coordinates, 714446, 684458, 689048, 17, nmax);

        simdtrf::compute_hrr_dl(buffer, coordinates, 719036, 686753, 691853, 17, nmax);

        simdtrf::compute_hrr_dm(buffer, coordinates, 723626, 689048, 694658, 17, nmax);

        simdtrf::compute_hrr_dm(buffer, coordinates, 729236, 691853, 698024, 17, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 734846, 701390, 707102, 17, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 739606, 704246, 710774, 17, nmax);

        simdtrf::compute_hrr_fk(buffer, coordinates, 744366, 707102, 714446, 17, nmax);

        simdtrf::compute_hrr_fk(buffer, coordinates, 750486, 710774, 719036, 17, nmax);

        simdtrf::compute_hrr_fl(buffer, coordinates, 756606, 714446, 723626, 17, nmax);

        simdtrf::compute_hrr_fl(buffer, coordinates, 764256, 719036, 729236, 17, nmax);

        simdtrf::compute_hrr_gi(buffer, coordinates, 771906, 734846, 744366, 17, nmax);

        simdtrf::compute_hrr_gi(buffer, coordinates, 779046, 739606, 750486, 17, nmax);

        simdtrf::compute_hrr_gk(buffer, coordinates, 786186, 744366, 756606, 17, nmax);

        simdtrf::compute_hrr_gk(buffer, coordinates, 795366, 750486, 764256, 17, nmax);

        simdtrf::compute_hrr_hi(buffer, coordinates, 804546, 771906, 786186, 17, nmax);

        simdtrf::compute_hrr_hi(buffer, coordinates, 814542, 779046, 795366, 17, nmax);

        simdtrf::transform_i_inner(buffer, 824538, 814542, 21, 17, nmax);

        simdtrf::transform_h_outer(values + n * npairs, nvalues, buffer, 824538, 221, nmax);

        simdtrf::transform_i_inner(buffer, 824538, 804546, 21, 17, nmax);

        simdtrf::transform_h_outer(values + 2431 * nvalues + n * npairs, nvalues, buffer, 824538,
                                   221, nmax);
    }

    for (size_t m = 0; m < 4862; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
