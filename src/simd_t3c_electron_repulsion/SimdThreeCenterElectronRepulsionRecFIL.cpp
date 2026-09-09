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


#include "SimdThreeCenterElectronRepulsionRecFIL.hpp"

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
#include "SimdTransferFI.hpp"
#include "SimdTransferPI.hpp"
#include "SimdTransferPK.hpp"
#include "SimdTransferPL.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformI.hpp"
#include "SimdTransformL.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_fil_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_fil_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 204407, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1547 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 204407, 175182, 9233, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 6, 3, 17,
                                                             ncols, fj, mu, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 25, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 28, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 31, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 34, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 37, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 40, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 43, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 46, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 49, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 52, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 55, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 58, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 61, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 64, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 67, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 70, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 73, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 76, 0, 3, 7, 8,
                                                                       25, 28, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 82, 0, 3, 8, 9,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 88, 0, 3, 9, 10,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 94, 0, 3, 10, 11,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 100, 0, 3, 11, 12,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 106, 0, 3, 12, 13,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 112, 0, 3, 13, 14,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 118, 0, 3, 14, 15,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 124, 0, 3, 15, 16,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 130, 0, 3, 16, 17,
                                                                       52, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 136, 0, 3, 17, 18,
                                                                       55, 58, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 142, 0, 3, 18, 19,
                                                                       58, 61, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 148, 0, 3, 19, 20,
                                                                       61, 64, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 154, 0, 3, 20, 21,
                                                                       64, 67, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 160, 0, 3, 21, 22,
                                                                       67, 70, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 166, 0, 3, 22, 23,
                                                                       70, 73, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 172, 0, 3, 25, 28,
                                                                       76, 82, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 182, 0, 3, 28, 31,
                                                                       82, 88, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 192, 0, 3, 31, 34,
                                                                       88, 94, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 202, 0, 3, 34, 37,
                                                                       94, 100, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 212, 0, 3, 37, 40,
                                                                       100, 106, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 222, 0, 3, 40, 43,
                                                                       106, 112, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 232, 0, 3, 43, 46,
                                                                       112, 118, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 242, 0, 3, 46, 49,
                                                                       118, 124, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 252, 0, 3, 49, 52,
                                                                       124, 130, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 262, 0, 3, 52, 55,
                                                                       130, 136, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 272, 0, 3, 55, 58,
                                                                       136, 142, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 282, 0, 3, 58, 61,
                                                                       142, 148, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 292, 0, 3, 61, 64,
                                                                       148, 154, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 302, 0, 3, 64, 67,
                                                                       154, 160, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 312, 0, 3, 67, 70,
                                                                       160, 166, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 322, 0, 3, 76, 82,
                                                                       172, 182, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 337, 0, 3, 82, 88,
                                                                       182, 192, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 352, 0, 3, 88, 94,
                                                                       192, 202, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 367, 0, 3, 94,
                                                                       100, 202, 212, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 382, 0, 3, 100,
                                                                       106, 212, 222, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 397, 0, 3, 106,
                                                                       112, 222, 232, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 412, 0, 3, 112,
                                                                       118, 232, 242, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 427, 0, 3, 118,
                                                                       124, 242, 252, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 442, 0, 3, 124,
                                                                       130, 252, 262, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 457, 0, 3, 130,
                                                                       136, 262, 272, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 472, 0, 3, 136,
                                                                       142, 272, 282, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 487, 0, 3, 142,
                                                                       148, 282, 292, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 502, 0, 3, 148,
                                                                       154, 292, 302, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 517, 0, 3, 154,
                                                                       160, 302, 312, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 532, 0, 3, 172,
                                                                       182, 322, 337, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 553, 0, 3, 182,
                                                                       192, 337, 352, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 574, 0, 3, 192,
                                                                       202, 352, 367, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 595, 0, 3, 202,
                                                                       212, 367, 382, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 616, 0, 3, 212,
                                                                       222, 382, 397, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 637, 0, 3, 222,
                                                                       232, 397, 412, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 658, 0, 3, 232,
                                                                       242, 412, 427, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 679, 0, 3, 242,
                                                                       252, 427, 442, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 700, 0, 3, 252,
                                                                       262, 442, 457, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 721, 0, 3, 262,
                                                                       272, 457, 472, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 742, 0, 3, 272,
                                                                       282, 472, 487, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 763, 0, 3, 282,
                                                                       292, 487, 502, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 784, 0, 3, 292,
                                                                       302, 502, 517, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 805, 0, 3, 322,
                                                                       337, 532, 553, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 833, 0, 3, 337,
                                                                       352, 553, 574, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 861, 0, 3, 352,
                                                                       367, 574, 595, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 889, 0, 3, 367,
                                                                       382, 595, 616, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 917, 0, 3, 382,
                                                                       397, 616, 637, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 945, 0, 3, 397,
                                                                       412, 637, 658, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 973, 0, 3, 412,
                                                                       427, 658, 679, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1001, 0, 3, 427,
                                                                       442, 679, 700, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1029, 0, 3, 442,
                                                                       457, 700, 721, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1057, 0, 3, 457,
                                                                       472, 721, 742, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1085, 0, 3, 472,
                                                                       487, 742, 763, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1113, 0, 3, 487,
                                                                       502, 763, 784, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1141, 0, 3, 532,
                                                                       553, 805, 833, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1177, 0, 3, 553,
                                                                       574, 833, 861, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1213, 0, 3, 574,
                                                                       595, 861, 889, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1249, 0, 3, 595,
                                                                       616, 889, 917, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1285, 0, 3, 616,
                                                                       637, 917, 945, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1321, 0, 3, 637,
                                                                       658, 945, 973, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1357, 0, 3, 658,
                                                                       679, 973, 1001, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1393, 0, 3, 679,
                                                                       700, 1001, 1029, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1429, 0, 3, 700,
                                                                       721, 1029, 1057, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1465, 0, 3, 721,
                                                                       742, 1057, 1085, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1501, 0, 3, 742,
                                                                       763, 1085, 1113, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1537, 0, 3, 805,
                                                                       833, 1141, 1177, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1582, 0, 3, 833,
                                                                       861, 1177, 1213, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1627, 0, 3, 861,
                                                                       889, 1213, 1249, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1672, 0, 3, 889,
                                                                       917, 1249, 1285, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1717, 0, 3, 917,
                                                                       945, 1285, 1321, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1762, 0, 3, 945,
                                                                       973, 1321, 1357, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1807, 0, 3, 973,
                                                                       1001, 1357, 1393, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1852, 0, 3, 1001,
                                                                       1029, 1393, 1429, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1897, 0, 3, 1029,
                                                                       1057, 1429, 1465, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1942, 0, 3, 1057,
                                                                       1085, 1465, 1501, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1987, 0, 3, 1141,
                                                                       1177, 1537, 1582, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2042, 0, 3, 1177,
                                                                       1213, 1582, 1627, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2097, 0, 3, 1213,
                                                                       1249, 1627, 1672, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2152, 0, 3, 1249,
                                                                       1285, 1672, 1717, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2207, 0, 3, 1285,
                                                                       1321, 1717, 1762, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2262, 0, 3, 1321,
                                                                       1357, 1762, 1807, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2317, 0, 3, 1357,
                                                                       1393, 1807, 1852, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2372, 0, 3, 1393,
                                                                       1429, 1852, 1897, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2427, 0, 3, 1429,
                                                                       1465, 1897, 1942, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2482, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2485, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2488, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2491, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2494, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2497, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2500, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2503, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2506, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2509, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2512, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2515, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2518, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2521, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2524, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2527, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2530, 3, 9, 31,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2539, 3, 10, 34,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2548, 3, 11, 37,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2557, 3, 12, 40,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2566, 3, 13, 43,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2575, 3, 14, 46,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2584, 3, 15, 49,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2593, 3, 16, 52,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2602, 3, 17, 55,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2611, 3, 18, 58,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2620, 3, 19, 61,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2629, 3, 20, 64,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2638, 3, 21, 67,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2647, 3, 22, 70,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2656, 3, 23, 73,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2665, 3, 31, 88,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2683, 3, 34, 94,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2701, 3, 37, 100,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2719, 3, 40, 106,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2737, 3, 43, 112,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2755, 3, 46, 118,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2773, 3, 49, 124,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2791, 3, 52, 130,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2809, 3, 55, 136,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2827, 3, 58, 142,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2845, 3, 61, 148,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2863, 3, 64, 154,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2881, 3, 67, 160,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2899, 3, 70, 166,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2917, 3, 88, 192,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2947, 3, 94, 202,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2977, 3, 100, 212,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3007, 3, 106, 222,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3037, 3, 112, 232,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3067, 3, 118, 242,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3097, 3, 124, 252,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3127, 3, 130, 262,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3157, 3, 136, 272,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3187, 3, 142, 282,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3217, 3, 148, 292,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3247, 3, 154, 302,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3277, 3, 160, 312,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3307, 3, 192, 352,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3352, 3, 202, 367,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3397, 3, 212, 382,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3442, 3, 222, 397,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3487, 3, 232, 412,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3532, 3, 242, 427,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3577, 3, 252, 442,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3622, 3, 262, 457,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3667, 3, 272, 472,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3712, 3, 282, 487,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3757, 3, 292, 502,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3802, 3, 302, 517,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3847, 3, 352, 574,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3910, 3, 367, 595,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3973, 3, 382, 616,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4036, 3, 397, 637,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4099, 3, 412, 658,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4162, 3, 427, 679,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4225, 3, 442, 700,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4288, 3, 457, 721,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4351, 3, 472, 742,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4414, 3, 487, 763,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4477, 3, 502, 784,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4540, 3, 574, 861,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4624, 3, 595, 889,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4708, 3, 616, 917,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4792, 3, 637, 945,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4876, 3, 658, 973,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4960, 3, 679,
                                                                       1001, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5044, 3, 700,
                                                                       1029, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5128, 3, 721,
                                                                       1057, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5212, 3, 742,
                                                                       1085, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5296, 3, 763,
                                                                       1113, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5380, 3, 861,
                                                                       1213, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5488, 3, 889,
                                                                       1249, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5596, 3, 917,
                                                                       1285, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5704, 3, 945,
                                                                       1321, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5812, 3, 973,
                                                                       1357, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5920, 3, 1001,
                                                                       1393, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6028, 3, 1029,
                                                                       1429, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6136, 3, 1057,
                                                                       1465, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6244, 3, 1085,
                                                                       1501, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6352, 3, 1213,
                                                                       1627, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6487, 3, 1249,
                                                                       1672, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6622, 3, 1285,
                                                                       1717, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6757, 3, 1321,
                                                                       1762, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6892, 3, 1357,
                                                                       1807, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 7027, 3, 1393,
                                                                       1852, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 7162, 3, 1429,
                                                                       1897, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 7297, 3, 1465,
                                                                       1942, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 7432, 3, 1627,
                                                                       2097, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 7597, 3, 1672,
                                                                       2152, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 7762, 3, 1717,
                                                                       2207, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 7927, 3, 1762,
                                                                       2262, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 8092, 3, 1807,
                                                                       2317, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 8257, 3, 1852,
                                                                       2372, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 8422, 3, 1897,
                                                                       2427, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8587, 3, 7, 8,
                                                                       2482, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8593, 3, 8, 9,
                                                                       2485, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8599, 3, 9, 10,
                                                                       2488, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8605, 3, 10, 11,
                                                                       2491, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8611, 3, 11, 12,
                                                                       2494, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8617, 3, 12, 13,
                                                                       2497, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8623, 3, 13, 14,
                                                                       2500, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8629, 3, 14, 15,
                                                                       2503, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8635, 3, 15, 16,
                                                                       2506, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8641, 3, 16, 17,
                                                                       2509, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8647, 3, 17, 18,
                                                                       2512, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8653, 3, 18, 19,
                                                                       2515, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8659, 3, 19, 20,
                                                                       2518, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8665, 3, 20, 21,
                                                                       2521, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8671, 3, 21, 22,
                                                                       2524, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8677, 3, 22, 23,
                                                                       2527, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8683, 0, 3, 8587,
                                                                       2482, 8593, 2530, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8701, 0, 3, 8593,
                                                                       2485, 8599, 2539, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8719, 0, 3, 8599,
                                                                       2488, 8605, 2548, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8737, 0, 3, 8605,
                                                                       2491, 8611, 2557, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8755, 0, 3, 8611,
                                                                       2494, 8617, 2566, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8773, 0, 3, 8617,
                                                                       2497, 8623, 2575, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8791, 0, 3, 8623,
                                                                       2500, 8629, 2584, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8809, 0, 3, 8629,
                                                                       2503, 8635, 2593, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8827, 0, 3, 8635,
                                                                       2506, 8641, 2602, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8845, 0, 3, 8641,
                                                                       2509, 8647, 2611, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8863, 0, 3, 8647,
                                                                       2512, 8653, 2620, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8881, 0, 3, 8653,
                                                                       2515, 8659, 2629, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8899, 0, 3, 8659,
                                                                       2518, 8665, 2638, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8917, 0, 3, 8665,
                                                                       2521, 8671, 2647, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8935, 0, 3, 8671,
                                                                       2524, 8677, 2656, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8953, 0, 3, 8683,
                                                                       2530, 8701, 76, 82, 2665,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8989, 0, 3, 8701,
                                                                       2539, 8719, 82, 88, 2683,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9025, 0, 3, 8719,
                                                                       2548, 8737, 88, 94, 2701,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9061, 0, 3, 8737,
                                                                       2557, 8755, 94, 100, 2719,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9097, 0, 3, 8755,
                                                                       2566, 8773, 100, 106,
                                                                       2737, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9133, 0, 3, 8773,
                                                                       2575, 8791, 106, 112,
                                                                       2755, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9169, 0, 3, 8791,
                                                                       2584, 8809, 112, 118,
                                                                       2773, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9205, 0, 3, 8809,
                                                                       2593, 8827, 118, 124,
                                                                       2791, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9241, 0, 3, 8827,
                                                                       2602, 8845, 124, 130,
                                                                       2809, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9277, 0, 3, 8845,
                                                                       2611, 8863, 130, 136,
                                                                       2827, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9313, 0, 3, 8863,
                                                                       2620, 8881, 136, 142,
                                                                       2845, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9349, 0, 3, 8881,
                                                                       2629, 8899, 142, 148,
                                                                       2863, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9385, 0, 3, 8899,
                                                                       2638, 8917, 148, 154,
                                                                       2881, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9421, 0, 3, 8917,
                                                                       2647, 8935, 154, 160,
                                                                       2899, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9457, 0, 3, 8953,
                                                                       2665, 8989, 172, 182,
                                                                       2917, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9517, 0, 3, 8989,
                                                                       2683, 9025, 182, 192,
                                                                       2947, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9577, 0, 3, 9025,
                                                                       2701, 9061, 192, 202,
                                                                       2977, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9637, 0, 3, 9061,
                                                                       2719, 9097, 202, 212,
                                                                       3007, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9697, 0, 3, 9097,
                                                                       2737, 9133, 212, 222,
                                                                       3037, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9757, 0, 3, 9133,
                                                                       2755, 9169, 222, 232,
                                                                       3067, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9817, 0, 3, 9169,
                                                                       2773, 9205, 232, 242,
                                                                       3097, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9877, 0, 3, 9205,
                                                                       2791, 9241, 242, 252,
                                                                       3127, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9937, 0, 3, 9241,
                                                                       2809, 9277, 252, 262,
                                                                       3157, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9997, 0, 3, 9277,
                                                                       2827, 9313, 262, 272,
                                                                       3187, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10057, 0, 3, 9313,
                                                                       2845, 9349, 272, 282,
                                                                       3217, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10117, 0, 3, 9349,
                                                                       2863, 9385, 282, 292,
                                                                       3247, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10177, 0, 3, 9385,
                                                                       2881, 9421, 292, 302,
                                                                       3277, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10237, 0, 3, 9457,
                                                                       2917, 9517, 322, 337,
                                                                       3307, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10327, 0, 3, 9517,
                                                                       2947, 9577, 337, 352,
                                                                       3352, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10417, 0, 3, 9577,
                                                                       2977, 9637, 352, 367,
                                                                       3397, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10507, 0, 3, 9637,
                                                                       3007, 9697, 367, 382,
                                                                       3442, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10597, 0, 3, 9697,
                                                                       3037, 9757, 382, 397,
                                                                       3487, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10687, 0, 3, 9757,
                                                                       3067, 9817, 397, 412,
                                                                       3532, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10777, 0, 3, 9817,
                                                                       3097, 9877, 412, 427,
                                                                       3577, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10867, 0, 3, 9877,
                                                                       3127, 9937, 427, 442,
                                                                       3622, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10957, 0, 3, 9937,
                                                                       3157, 9997, 442, 457,
                                                                       3667, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11047, 0, 3, 9997,
                                                                       3187, 10057, 457, 472,
                                                                       3712, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11137, 0, 3,
                                                                       10057, 3217, 10117, 472,
                                                                       487, 3757, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11227, 0, 3,
                                                                       10117, 3247, 10177, 487,
                                                                       502, 3802, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11317, 0, 3,
                                                                       10237, 3307, 10327, 532,
                                                                       553, 3847, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11443, 0, 3,
                                                                       10327, 3352, 10417, 553,
                                                                       574, 3910, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11569, 0, 3,
                                                                       10417, 3397, 10507, 574,
                                                                       595, 3973, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11695, 0, 3,
                                                                       10507, 3442, 10597, 595,
                                                                       616, 4036, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11821, 0, 3,
                                                                       10597, 3487, 10687, 616,
                                                                       637, 4099, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11947, 0, 3,
                                                                       10687, 3532, 10777, 637,
                                                                       658, 4162, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 12073, 0, 3,
                                                                       10777, 3577, 10867, 658,
                                                                       679, 4225, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 12199, 0, 3,
                                                                       10867, 3622, 10957, 679,
                                                                       700, 4288, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 12325, 0, 3,
                                                                       10957, 3667, 11047, 700,
                                                                       721, 4351, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 12451, 0, 3,
                                                                       11047, 3712, 11137, 721,
                                                                       742, 4414, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 12577, 0, 3,
                                                                       11137, 3757, 11227, 742,
                                                                       763, 4477, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 12703, 0, 3,
                                                                       11317, 3847, 11443, 805,
                                                                       833, 4540, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 12871, 0, 3,
                                                                       11443, 3910, 11569, 833,
                                                                       861, 4624, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 13039, 0, 3,
                                                                       11569, 3973, 11695, 861,
                                                                       889, 4708, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 13207, 0, 3,
                                                                       11695, 4036, 11821, 889,
                                                                       917, 4792, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 13375, 0, 3,
                                                                       11821, 4099, 11947, 917,
                                                                       945, 4876, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 13543, 0, 3,
                                                                       11947, 4162, 12073, 945,
                                                                       973, 4960, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 13711, 0, 3,
                                                                       12073, 4225, 12199, 973,
                                                                       1001, 5044, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 13879, 0, 3,
                                                                       12199, 4288, 12325, 1001,
                                                                       1029, 5128, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 14047, 0, 3,
                                                                       12325, 4351, 12451, 1029,
                                                                       1057, 5212, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 14215, 0, 3,
                                                                       12451, 4414, 12577, 1057,
                                                                       1085, 5296, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 14383, 0, 3,
                                                                       12703, 4540, 12871, 1141,
                                                                       1177, 5380, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 14599, 0, 3,
                                                                       12871, 4624, 13039, 1177,
                                                                       1213, 5488, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 14815, 0, 3,
                                                                       13039, 4708, 13207, 1213,
                                                                       1249, 5596, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 15031, 0, 3,
                                                                       13207, 4792, 13375, 1249,
                                                                       1285, 5704, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 15247, 0, 3,
                                                                       13375, 4876, 13543, 1285,
                                                                       1321, 5812, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 15463, 0, 3,
                                                                       13543, 4960, 13711, 1321,
                                                                       1357, 5920, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 15679, 0, 3,
                                                                       13711, 5044, 13879, 1357,
                                                                       1393, 6028, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 15895, 0, 3,
                                                                       13879, 5128, 14047, 1393,
                                                                       1429, 6136, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 16111, 0, 3,
                                                                       14047, 5212, 14215, 1429,
                                                                       1465, 6244, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 16327, 0, 3,
                                                                       14383, 5380, 14599, 1537,
                                                                       1582, 6352, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 16597, 0, 3,
                                                                       14599, 5488, 14815, 1582,
                                                                       1627, 6487, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 16867, 0, 3,
                                                                       14815, 5596, 15031, 1627,
                                                                       1672, 6622, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 17137, 0, 3,
                                                                       15031, 5704, 15247, 1672,
                                                                       1717, 6757, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 17407, 0, 3,
                                                                       15247, 5812, 15463, 1717,
                                                                       1762, 6892, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 17677, 0, 3,
                                                                       15463, 5920, 15679, 1762,
                                                                       1807, 7027, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 17947, 0, 3,
                                                                       15679, 6028, 15895, 1807,
                                                                       1852, 7162, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 18217, 0, 3,
                                                                       15895, 6136, 16111, 1852,
                                                                       1897, 7297, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 18487, 0, 3,
                                                                       16327, 6352, 16597, 1987,
                                                                       2042, 7432, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 18817, 0, 3,
                                                                       16597, 6487, 16867, 2042,
                                                                       2097, 7597, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 19147, 0, 3,
                                                                       16867, 6622, 17137, 2097,
                                                                       2152, 7762, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 19477, 0, 3,
                                                                       17137, 6757, 17407, 2152,
                                                                       2207, 7927, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 19807, 0, 3,
                                                                       17407, 6892, 17677, 2207,
                                                                       2262, 8092, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 20137, 0, 3,
                                                                       17677, 7027, 17947, 2262,
                                                                       2317, 8257, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 20467, 0, 3,
                                                                       17947, 7162, 18217, 2317,
                                                                       2372, 8422, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 20797, 3, 2482,
                                                                       2485, 8599, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 20807, 3, 2485,
                                                                       2488, 8605, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 20817, 3, 2488,
                                                                       2491, 8611, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 20827, 3, 2491,
                                                                       2494, 8617, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 20837, 3, 2494,
                                                                       2497, 8623, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 20847, 3, 2497,
                                                                       2500, 8629, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 20857, 3, 2500,
                                                                       2503, 8635, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 20867, 3, 2503,
                                                                       2506, 8641, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 20877, 3, 2506,
                                                                       2509, 8647, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 20887, 3, 2509,
                                                                       2512, 8653, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 20897, 3, 2512,
                                                                       2515, 8659, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 20907, 3, 2515,
                                                                       2518, 8665, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 20917, 3, 2518,
                                                                       2521, 8671, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 20927, 3, 2521,
                                                                       2524, 8677, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 20937, 0, 3,
                                                                       20797, 8599, 20807, 8719,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 20967, 0, 3,
                                                                       20807, 8605, 20817, 8737,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 20997, 0, 3,
                                                                       20817, 8611, 20827, 8755,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 21027, 0, 3,
                                                                       20827, 8617, 20837, 8773,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 21057, 0, 3,
                                                                       20837, 8623, 20847, 8791,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 21087, 0, 3,
                                                                       20847, 8629, 20857, 8809,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 21117, 0, 3,
                                                                       20857, 8635, 20867, 8827,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 21147, 0, 3,
                                                                       20867, 8641, 20877, 8845,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 21177, 0, 3,
                                                                       20877, 8647, 20887, 8863,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 21207, 0, 3,
                                                                       20887, 8653, 20897, 8881,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 21237, 0, 3,
                                                                       20897, 8659, 20907, 8899,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 21267, 0, 3,
                                                                       20907, 8665, 20917, 8917,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 21297, 0, 3,
                                                                       20917, 8671, 20927, 8935,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 21327, 0, 3,
                                                                       20937, 8719, 20967, 2665,
                                                                       2683, 9025, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 21387, 0, 3,
                                                                       20967, 8737, 20997, 2683,
                                                                       2701, 9061, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 21447, 0, 3,
                                                                       20997, 8755, 21027, 2701,
                                                                       2719, 9097, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 21507, 0, 3,
                                                                       21027, 8773, 21057, 2719,
                                                                       2737, 9133, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 21567, 0, 3,
                                                                       21057, 8791, 21087, 2737,
                                                                       2755, 9169, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 21627, 0, 3,
                                                                       21087, 8809, 21117, 2755,
                                                                       2773, 9205, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 21687, 0, 3,
                                                                       21117, 8827, 21147, 2773,
                                                                       2791, 9241, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 21747, 0, 3,
                                                                       21147, 8845, 21177, 2791,
                                                                       2809, 9277, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 21807, 0, 3,
                                                                       21177, 8863, 21207, 2809,
                                                                       2827, 9313, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 21867, 0, 3,
                                                                       21207, 8881, 21237, 2827,
                                                                       2845, 9349, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 21927, 0, 3,
                                                                       21237, 8899, 21267, 2845,
                                                                       2863, 9385, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 21987, 0, 3,
                                                                       21267, 8917, 21297, 2863,
                                                                       2881, 9421, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 22047, 0, 3,
                                                                       21327, 9025, 21387, 2917,
                                                                       2947, 9577, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 22147, 0, 3,
                                                                       21387, 9061, 21447, 2947,
                                                                       2977, 9637, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 22247, 0, 3,
                                                                       21447, 9097, 21507, 2977,
                                                                       3007, 9697, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 22347, 0, 3,
                                                                       21507, 9133, 21567, 3007,
                                                                       3037, 9757, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 22447, 0, 3,
                                                                       21567, 9169, 21627, 3037,
                                                                       3067, 9817, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 22547, 0, 3,
                                                                       21627, 9205, 21687, 3067,
                                                                       3097, 9877, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 22647, 0, 3,
                                                                       21687, 9241, 21747, 3097,
                                                                       3127, 9937, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 22747, 0, 3,
                                                                       21747, 9277, 21807, 3127,
                                                                       3157, 9997, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 22847, 0, 3,
                                                                       21807, 9313, 21867, 3157,
                                                                       3187, 10057, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 22947, 0, 3,
                                                                       21867, 9349, 21927, 3187,
                                                                       3217, 10117, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 23047, 0, 3,
                                                                       21927, 9385, 21987, 3217,
                                                                       3247, 10177, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 23147, 0, 3,
                                                                       22047, 9577, 22147, 3307,
                                                                       3352, 10417, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 23297, 0, 3,
                                                                       22147, 9637, 22247, 3352,
                                                                       3397, 10507, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 23447, 0, 3,
                                                                       22247, 9697, 22347, 3397,
                                                                       3442, 10597, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 23597, 0, 3,
                                                                       22347, 9757, 22447, 3442,
                                                                       3487, 10687, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 23747, 0, 3,
                                                                       22447, 9817, 22547, 3487,
                                                                       3532, 10777, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 23897, 0, 3,
                                                                       22547, 9877, 22647, 3532,
                                                                       3577, 10867, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 24047, 0, 3,
                                                                       22647, 9937, 22747, 3577,
                                                                       3622, 10957, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 24197, 0, 3,
                                                                       22747, 9997, 22847, 3622,
                                                                       3667, 11047, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 24347, 0, 3,
                                                                       22847, 10057, 22947, 3667,
                                                                       3712, 11137, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 24497, 0, 3,
                                                                       22947, 10117, 23047, 3712,
                                                                       3757, 11227, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 24647, 0, 3,
                                                                       23147, 10417, 23297, 3847,
                                                                       3910, 11569, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 24857, 0, 3,
                                                                       23297, 10507, 23447, 3910,
                                                                       3973, 11695, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 25067, 0, 3,
                                                                       23447, 10597, 23597, 3973,
                                                                       4036, 11821, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 25277, 0, 3,
                                                                       23597, 10687, 23747, 4036,
                                                                       4099, 11947, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 25487, 0, 3,
                                                                       23747, 10777, 23897, 4099,
                                                                       4162, 12073, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 25697, 0, 3,
                                                                       23897, 10867, 24047, 4162,
                                                                       4225, 12199, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 25907, 0, 3,
                                                                       24047, 10957, 24197, 4225,
                                                                       4288, 12325, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 26117, 0, 3,
                                                                       24197, 11047, 24347, 4288,
                                                                       4351, 12451, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 26327, 0, 3,
                                                                       24347, 11137, 24497, 4351,
                                                                       4414, 12577, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 26537, 0, 3,
                                                                       24647, 11569, 24857, 4540,
                                                                       4624, 13039, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 26817, 0, 3,
                                                                       24857, 11695, 25067, 4624,
                                                                       4708, 13207, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 27097, 0, 3,
                                                                       25067, 11821, 25277, 4708,
                                                                       4792, 13375, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 27377, 0, 3,
                                                                       25277, 11947, 25487, 4792,
                                                                       4876, 13543, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 27657, 0, 3,
                                                                       25487, 12073, 25697, 4876,
                                                                       4960, 13711, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 27937, 0, 3,
                                                                       25697, 12199, 25907, 4960,
                                                                       5044, 13879, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 28217, 0, 3,
                                                                       25907, 12325, 26117, 5044,
                                                                       5128, 14047, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 28497, 0, 3,
                                                                       26117, 12451, 26327, 5128,
                                                                       5212, 14215, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 28777, 0, 3,
                                                                       26537, 13039, 26817, 5380,
                                                                       5488, 14815, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 29137, 0, 3,
                                                                       26817, 13207, 27097, 5488,
                                                                       5596, 15031, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 29497, 0, 3,
                                                                       27097, 13375, 27377, 5596,
                                                                       5704, 15247, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 29857, 0, 3,
                                                                       27377, 13543, 27657, 5704,
                                                                       5812, 15463, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 30217, 0, 3,
                                                                       27657, 13711, 27937, 5812,
                                                                       5920, 15679, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 30577, 0, 3,
                                                                       27937, 13879, 28217, 5920,
                                                                       6028, 15895, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 30937, 0, 3,
                                                                       28217, 14047, 28497, 6028,
                                                                       6136, 16111, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 31297, 0, 3,
                                                                       28777, 14815, 29137, 6352,
                                                                       6487, 16867, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 31747, 0, 3,
                                                                       29137, 15031, 29497, 6487,
                                                                       6622, 17137, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 32197, 0, 3,
                                                                       29497, 15247, 29857, 6622,
                                                                       6757, 17407, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 32647, 0, 3,
                                                                       29857, 15463, 30217, 6757,
                                                                       6892, 17677, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 33097, 0, 3,
                                                                       30217, 15679, 30577, 6892,
                                                                       7027, 17947, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 33547, 0, 3,
                                                                       30577, 15895, 30937, 7027,
                                                                       7162, 18217, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 33997, 0, 3,
                                                                       31297, 16867, 31747, 7432,
                                                                       7597, 19147, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 34547, 0, 3,
                                                                       31747, 17137, 32197, 7597,
                                                                       7762, 19477, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 35097, 0, 3,
                                                                       32197, 17407, 32647, 7762,
                                                                       7927, 19807, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 35647, 0, 3,
                                                                       32647, 17677, 33097, 7927,
                                                                       8092, 20137, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 36197, 0, 3,
                                                                       33097, 17947, 33547, 8092,
                                                                       8257, 20467, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 36747, 3, 8587,
                                                                       8593, 20797, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 36762, 3, 8593,
                                                                       8599, 20807, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 36777, 3, 8599,
                                                                       8605, 20817, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 36792, 3, 8605,
                                                                       8611, 20827, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 36807, 3, 8611,
                                                                       8617, 20837, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 36822, 3, 8617,
                                                                       8623, 20847, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 36837, 3, 8623,
                                                                       8629, 20857, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 36852, 3, 8629,
                                                                       8635, 20867, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 36867, 3, 8635,
                                                                       8641, 20877, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 36882, 3, 8641,
                                                                       8647, 20887, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 36897, 3, 8647,
                                                                       8653, 20897, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 36912, 3, 8653,
                                                                       8659, 20907, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 36927, 3, 8659,
                                                                       8665, 20917, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 36942, 3, 8665,
                                                                       8671, 20927, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 36957, 0, 3,
                                                                       36747, 20797, 36762, 8683,
                                                                       8701, 20937, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 37002, 0, 3,
                                                                       36762, 20807, 36777, 8701,
                                                                       8719, 20967, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 37047, 0, 3,
                                                                       36777, 20817, 36792, 8719,
                                                                       8737, 20997, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 37092, 0, 3,
                                                                       36792, 20827, 36807, 8737,
                                                                       8755, 21027, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 37137, 0, 3,
                                                                       36807, 20837, 36822, 8755,
                                                                       8773, 21057, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 37182, 0, 3,
                                                                       36822, 20847, 36837, 8773,
                                                                       8791, 21087, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 37227, 0, 3,
                                                                       36837, 20857, 36852, 8791,
                                                                       8809, 21117, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 37272, 0, 3,
                                                                       36852, 20867, 36867, 8809,
                                                                       8827, 21147, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 37317, 0, 3,
                                                                       36867, 20877, 36882, 8827,
                                                                       8845, 21177, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 37362, 0, 3,
                                                                       36882, 20887, 36897, 8845,
                                                                       8863, 21207, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 37407, 0, 3,
                                                                       36897, 20897, 36912, 8863,
                                                                       8881, 21237, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 37452, 0, 3,
                                                                       36912, 20907, 36927, 8881,
                                                                       8899, 21267, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 37497, 0, 3,
                                                                       36927, 20917, 36942, 8899,
                                                                       8917, 21297, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 37542, 0, 3,
                                                                       36957, 20937, 37002, 8953,
                                                                       8989, 21327, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 37632, 0, 3,
                                                                       37002, 20967, 37047, 8989,
                                                                       9025, 21387, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 37722, 0, 3,
                                                                       37047, 20997, 37092, 9025,
                                                                       9061, 21447, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 37812, 0, 3,
                                                                       37092, 21027, 37137, 9061,
                                                                       9097, 21507, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 37902, 0, 3,
                                                                       37137, 21057, 37182, 9097,
                                                                       9133, 21567, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 37992, 0, 3,
                                                                       37182, 21087, 37227, 9133,
                                                                       9169, 21627, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 38082, 0, 3,
                                                                       37227, 21117, 37272, 9169,
                                                                       9205, 21687, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 38172, 0, 3,
                                                                       37272, 21147, 37317, 9205,
                                                                       9241, 21747, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 38262, 0, 3,
                                                                       37317, 21177, 37362, 9241,
                                                                       9277, 21807, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 38352, 0, 3,
                                                                       37362, 21207, 37407, 9277,
                                                                       9313, 21867, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 38442, 0, 3,
                                                                       37407, 21237, 37452, 9313,
                                                                       9349, 21927, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 38532, 0, 3,
                                                                       37452, 21267, 37497, 9349,
                                                                       9385, 21987, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 38622, 0, 3,
                                                                       37542, 21327, 37632, 9457,
                                                                       9517, 22047, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 38772, 0, 3,
                                                                       37632, 21387, 37722, 9517,
                                                                       9577, 22147, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 38922, 0, 3,
                                                                       37722, 21447, 37812, 9577,
                                                                       9637, 22247, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 39072, 0, 3,
                                                                       37812, 21507, 37902, 9637,
                                                                       9697, 22347, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 39222, 0, 3,
                                                                       37902, 21567, 37992, 9697,
                                                                       9757, 22447, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 39372, 0, 3,
                                                                       37992, 21627, 38082, 9757,
                                                                       9817, 22547, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 39522, 0, 3,
                                                                       38082, 21687, 38172, 9817,
                                                                       9877, 22647, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 39672, 0, 3,
                                                                       38172, 21747, 38262, 9877,
                                                                       9937, 22747, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 39822, 0, 3,
                                                                       38262, 21807, 38352, 9937,
                                                                       9997, 22847, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 39972, 0, 3,
                                                                       38352, 21867, 38442, 9997,
                                                                       10057, 22947, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 40122, 0, 3,
                                                                       38442, 21927, 38532,
                                                                       10057, 10117, 23047,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 40272, 0, 3,
                                                                       38622, 22047, 38772,
                                                                       10237, 10327, 23147,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 40497, 0, 3,
                                                                       38772, 22147, 38922,
                                                                       10327, 10417, 23297,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 40722, 0, 3,
                                                                       38922, 22247, 39072,
                                                                       10417, 10507, 23447,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 40947, 0, 3,
                                                                       39072, 22347, 39222,
                                                                       10507, 10597, 23597,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 41172, 0, 3,
                                                                       39222, 22447, 39372,
                                                                       10597, 10687, 23747,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 41397, 0, 3,
                                                                       39372, 22547, 39522,
                                                                       10687, 10777, 23897,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 41622, 0, 3,
                                                                       39522, 22647, 39672,
                                                                       10777, 10867, 24047,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 41847, 0, 3,
                                                                       39672, 22747, 39822,
                                                                       10867, 10957, 24197,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 42072, 0, 3,
                                                                       39822, 22847, 39972,
                                                                       10957, 11047, 24347,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 42297, 0, 3,
                                                                       39972, 22947, 40122,
                                                                       11047, 11137, 24497,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 42522, 0, 3,
                                                                       40272, 23147, 40497,
                                                                       11317, 11443, 24647,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 42837, 0, 3,
                                                                       40497, 23297, 40722,
                                                                       11443, 11569, 24857,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 43152, 0, 3,
                                                                       40722, 23447, 40947,
                                                                       11569, 11695, 25067,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 43467, 0, 3,
                                                                       40947, 23597, 41172,
                                                                       11695, 11821, 25277,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 43782, 0, 3,
                                                                       41172, 23747, 41397,
                                                                       11821, 11947, 25487,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 44097, 0, 3,
                                                                       41397, 23897, 41622,
                                                                       11947, 12073, 25697,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 44412, 0, 3,
                                                                       41622, 24047, 41847,
                                                                       12073, 12199, 25907,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 44727, 0, 3,
                                                                       41847, 24197, 42072,
                                                                       12199, 12325, 26117,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 45042, 0, 3,
                                                                       42072, 24347, 42297,
                                                                       12325, 12451, 26327,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 45357, 0, 3,
                                                                       42522, 24647, 42837,
                                                                       12703, 12871, 26537,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 45777, 0, 3,
                                                                       42837, 24857, 43152,
                                                                       12871, 13039, 26817,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 46197, 0, 3,
                                                                       43152, 25067, 43467,
                                                                       13039, 13207, 27097,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 46617, 0, 3,
                                                                       43467, 25277, 43782,
                                                                       13207, 13375, 27377,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 47037, 0, 3,
                                                                       43782, 25487, 44097,
                                                                       13375, 13543, 27657,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 47457, 0, 3,
                                                                       44097, 25697, 44412,
                                                                       13543, 13711, 27937,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 47877, 0, 3,
                                                                       44412, 25907, 44727,
                                                                       13711, 13879, 28217,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 48297, 0, 3,
                                                                       44727, 26117, 45042,
                                                                       13879, 14047, 28497,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 48717, 0, 3,
                                                                       45357, 26537, 45777,
                                                                       14383, 14599, 28777,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 49257, 0, 3,
                                                                       45777, 26817, 46197,
                                                                       14599, 14815, 29137,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 49797, 0, 3,
                                                                       46197, 27097, 46617,
                                                                       14815, 15031, 29497,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 50337, 0, 3,
                                                                       46617, 27377, 47037,
                                                                       15031, 15247, 29857,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 50877, 0, 3,
                                                                       47037, 27657, 47457,
                                                                       15247, 15463, 30217,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 51417, 0, 3,
                                                                       47457, 27937, 47877,
                                                                       15463, 15679, 30577,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 51957, 0, 3,
                                                                       47877, 28217, 48297,
                                                                       15679, 15895, 30937,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 52497, 0, 3,
                                                                       48717, 28777, 49257,
                                                                       16327, 16597, 31297,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 53172, 0, 3,
                                                                       49257, 29137, 49797,
                                                                       16597, 16867, 31747,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 53847, 0, 3,
                                                                       49797, 29497, 50337,
                                                                       16867, 17137, 32197,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 54522, 0, 3,
                                                                       50337, 29857, 50877,
                                                                       17137, 17407, 32647,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 55197, 0, 3,
                                                                       50877, 30217, 51417,
                                                                       17407, 17677, 33097,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 55872, 0, 3,
                                                                       51417, 30577, 51957,
                                                                       17677, 17947, 33547,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 56547, 0, 3,
                                                                       52497, 31297, 53172,
                                                                       18487, 18817, 33997,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 57372, 0, 3,
                                                                       53172, 31747, 53847,
                                                                       18817, 19147, 34547,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 58197, 0, 3,
                                                                       53847, 32197, 54522,
                                                                       19147, 19477, 35097,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 59022, 0, 3,
                                                                       54522, 32647, 55197,
                                                                       19477, 19807, 35647,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 59847, 0, 3,
                                                                       55197, 33097, 55872,
                                                                       19807, 20137, 36197,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 60672, 3, 20797,
                                                                       20807, 36777, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 60693, 3, 20807,
                                                                       20817, 36792, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 60714, 3, 20817,
                                                                       20827, 36807, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 60735, 3, 20827,
                                                                       20837, 36822, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 60756, 3, 20837,
                                                                       20847, 36837, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 60777, 3, 20847,
                                                                       20857, 36852, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 60798, 3, 20857,
                                                                       20867, 36867, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 60819, 3, 20867,
                                                                       20877, 36882, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 60840, 3, 20877,
                                                                       20887, 36897, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 60861, 3, 20887,
                                                                       20897, 36912, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 60882, 3, 20897,
                                                                       20907, 36927, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 60903, 3, 20907,
                                                                       20917, 36942, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 60924, 0, 3,
                                                                       60672, 36777, 60693,
                                                                       20937, 20967, 37047,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 60987, 0, 3,
                                                                       60693, 36792, 60714,
                                                                       20967, 20997, 37092,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 61050, 0, 3,
                                                                       60714, 36807, 60735,
                                                                       20997, 21027, 37137,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 61113, 0, 3,
                                                                       60735, 36822, 60756,
                                                                       21027, 21057, 37182,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 61176, 0, 3,
                                                                       60756, 36837, 60777,
                                                                       21057, 21087, 37227,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 61239, 0, 3,
                                                                       60777, 36852, 60798,
                                                                       21087, 21117, 37272,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 61302, 0, 3,
                                                                       60798, 36867, 60819,
                                                                       21117, 21147, 37317,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 61365, 0, 3,
                                                                       60819, 36882, 60840,
                                                                       21147, 21177, 37362,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 61428, 0, 3,
                                                                       60840, 36897, 60861,
                                                                       21177, 21207, 37407,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 61491, 0, 3,
                                                                       60861, 36912, 60882,
                                                                       21207, 21237, 37452,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 61554, 0, 3,
                                                                       60882, 36927, 60903,
                                                                       21237, 21267, 37497,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 61617, 0, 3,
                                                                       60924, 37047, 60987,
                                                                       21327, 21387, 37722,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 61743, 0, 3,
                                                                       60987, 37092, 61050,
                                                                       21387, 21447, 37812,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 61869, 0, 3,
                                                                       61050, 37137, 61113,
                                                                       21447, 21507, 37902,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 61995, 0, 3,
                                                                       61113, 37182, 61176,
                                                                       21507, 21567, 37992,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 62121, 0, 3,
                                                                       61176, 37227, 61239,
                                                                       21567, 21627, 38082,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 62247, 0, 3,
                                                                       61239, 37272, 61302,
                                                                       21627, 21687, 38172,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 62373, 0, 3,
                                                                       61302, 37317, 61365,
                                                                       21687, 21747, 38262,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 62499, 0, 3,
                                                                       61365, 37362, 61428,
                                                                       21747, 21807, 38352,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 62625, 0, 3,
                                                                       61428, 37407, 61491,
                                                                       21807, 21867, 38442,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 62751, 0, 3,
                                                                       61491, 37452, 61554,
                                                                       21867, 21927, 38532,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 62877, 0, 3,
                                                                       61617, 37722, 61743,
                                                                       22047, 22147, 38922,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 63087, 0, 3,
                                                                       61743, 37812, 61869,
                                                                       22147, 22247, 39072,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 63297, 0, 3,
                                                                       61869, 37902, 61995,
                                                                       22247, 22347, 39222,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 63507, 0, 3,
                                                                       61995, 37992, 62121,
                                                                       22347, 22447, 39372,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 63717, 0, 3,
                                                                       62121, 38082, 62247,
                                                                       22447, 22547, 39522,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 63927, 0, 3,
                                                                       62247, 38172, 62373,
                                                                       22547, 22647, 39672,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 64137, 0, 3,
                                                                       62373, 38262, 62499,
                                                                       22647, 22747, 39822,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 64347, 0, 3,
                                                                       62499, 38352, 62625,
                                                                       22747, 22847, 39972,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 64557, 0, 3,
                                                                       62625, 38442, 62751,
                                                                       22847, 22947, 40122,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 64767, 0, 3,
                                                                       62877, 38922, 63087,
                                                                       23147, 23297, 40722,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 65082, 0, 3,
                                                                       63087, 39072, 63297,
                                                                       23297, 23447, 40947,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 65397, 0, 3,
                                                                       63297, 39222, 63507,
                                                                       23447, 23597, 41172,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 65712, 0, 3,
                                                                       63507, 39372, 63717,
                                                                       23597, 23747, 41397,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 66027, 0, 3,
                                                                       63717, 39522, 63927,
                                                                       23747, 23897, 41622,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 66342, 0, 3,
                                                                       63927, 39672, 64137,
                                                                       23897, 24047, 41847,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 66657, 0, 3,
                                                                       64137, 39822, 64347,
                                                                       24047, 24197, 42072,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 66972, 0, 3,
                                                                       64347, 39972, 64557,
                                                                       24197, 24347, 42297,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 67287, 0, 3,
                                                                       64767, 40722, 65082,
                                                                       24647, 24857, 43152,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 67728, 0, 3,
                                                                       65082, 40947, 65397,
                                                                       24857, 25067, 43467,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 68169, 0, 3,
                                                                       65397, 41172, 65712,
                                                                       25067, 25277, 43782,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 68610, 0, 3,
                                                                       65712, 41397, 66027,
                                                                       25277, 25487, 44097,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 69051, 0, 3,
                                                                       66027, 41622, 66342,
                                                                       25487, 25697, 44412,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 69492, 0, 3,
                                                                       66342, 41847, 66657,
                                                                       25697, 25907, 44727,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 69933, 0, 3,
                                                                       66657, 42072, 66972,
                                                                       25907, 26117, 45042,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 70374, 0, 3,
                                                                       67287, 43152, 67728,
                                                                       26537, 26817, 46197,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 70962, 0, 3,
                                                                       67728, 43467, 68169,
                                                                       26817, 27097, 46617,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 71550, 0, 3,
                                                                       68169, 43782, 68610,
                                                                       27097, 27377, 47037,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 72138, 0, 3,
                                                                       68610, 44097, 69051,
                                                                       27377, 27657, 47457,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 72726, 0, 3,
                                                                       69051, 44412, 69492,
                                                                       27657, 27937, 47877,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 73314, 0, 3,
                                                                       69492, 44727, 69933,
                                                                       27937, 28217, 48297,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 73902, 0, 3,
                                                                       70374, 46197, 70962,
                                                                       28777, 29137, 49797,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 74658, 0, 3,
                                                                       70962, 46617, 71550,
                                                                       29137, 29497, 50337,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 75414, 0, 3,
                                                                       71550, 47037, 72138,
                                                                       29497, 29857, 50877,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 76170, 0, 3,
                                                                       72138, 47457, 72726,
                                                                       29857, 30217, 51417,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 76926, 0, 3,
                                                                       72726, 47877, 73314,
                                                                       30217, 30577, 51957,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 77682, 0, 3,
                                                                       73902, 49797, 74658,
                                                                       31297, 31747, 53847,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 78627, 0, 3,
                                                                       74658, 50337, 75414,
                                                                       31747, 32197, 54522,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 79572, 0, 3,
                                                                       75414, 50877, 76170,
                                                                       32197, 32647, 55197,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 80517, 0, 3,
                                                                       76170, 51417, 76926,
                                                                       32647, 33097, 55872,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 81462, 0, 3,
                                                                       77682, 53847, 78627,
                                                                       33997, 34547, 58197,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 82617, 0, 3,
                                                                       78627, 54522, 79572,
                                                                       34547, 35097, 59022,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 83772, 0, 3,
                                                                       79572, 55197, 80517,
                                                                       35097, 35647, 59847,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 84927, 3, 36747,
                                                                       36762, 60672, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 84955, 3, 36762,
                                                                       36777, 60693, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 84983, 3, 36777,
                                                                       36792, 60714, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 85011, 3, 36792,
                                                                       36807, 60735, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 85039, 3, 36807,
                                                                       36822, 60756, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 85067, 3, 36822,
                                                                       36837, 60777, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 85095, 3, 36837,
                                                                       36852, 60798, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 85123, 3, 36852,
                                                                       36867, 60819, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 85151, 3, 36867,
                                                                       36882, 60840, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 85179, 3, 36882,
                                                                       36897, 60861, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 85207, 3, 36897,
                                                                       36912, 60882, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 85235, 3, 36912,
                                                                       36927, 60903, ncols,
                                                                       gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 85263, 0, 3,
                                                                       84927, 60672, 84955,
                                                                       36957, 37002, 60924,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 85347, 0, 3,
                                                                       84955, 60693, 84983,
                                                                       37002, 37047, 60987,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 85431, 0, 3,
                                                                       84983, 60714, 85011,
                                                                       37047, 37092, 61050,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 85515, 0, 3,
                                                                       85011, 60735, 85039,
                                                                       37092, 37137, 61113,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 85599, 0, 3,
                                                                       85039, 60756, 85067,
                                                                       37137, 37182, 61176,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 85683, 0, 3,
                                                                       85067, 60777, 85095,
                                                                       37182, 37227, 61239,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 85767, 0, 3,
                                                                       85095, 60798, 85123,
                                                                       37227, 37272, 61302,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 85851, 0, 3,
                                                                       85123, 60819, 85151,
                                                                       37272, 37317, 61365,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 85935, 0, 3,
                                                                       85151, 60840, 85179,
                                                                       37317, 37362, 61428,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 86019, 0, 3,
                                                                       85179, 60861, 85207,
                                                                       37362, 37407, 61491,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 86103, 0, 3,
                                                                       85207, 60882, 85235,
                                                                       37407, 37452, 61554,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 86187, 0, 3,
                                                                       85263, 60924, 85347,
                                                                       37542, 37632, 61617,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 86355, 0, 3,
                                                                       85347, 60987, 85431,
                                                                       37632, 37722, 61743,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 86523, 0, 3,
                                                                       85431, 61050, 85515,
                                                                       37722, 37812, 61869,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 86691, 0, 3,
                                                                       85515, 61113, 85599,
                                                                       37812, 37902, 61995,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 86859, 0, 3,
                                                                       85599, 61176, 85683,
                                                                       37902, 37992, 62121,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 87027, 0, 3,
                                                                       85683, 61239, 85767,
                                                                       37992, 38082, 62247,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 87195, 0, 3,
                                                                       85767, 61302, 85851,
                                                                       38082, 38172, 62373,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 87363, 0, 3,
                                                                       85851, 61365, 85935,
                                                                       38172, 38262, 62499,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 87531, 0, 3,
                                                                       85935, 61428, 86019,
                                                                       38262, 38352, 62625,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 87699, 0, 3,
                                                                       86019, 61491, 86103,
                                                                       38352, 38442, 62751,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 87867, 0, 3,
                                                                       86187, 61617, 86355,
                                                                       38622, 38772, 62877,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 88147, 0, 3,
                                                                       86355, 61743, 86523,
                                                                       38772, 38922, 63087,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 88427, 0, 3,
                                                                       86523, 61869, 86691,
                                                                       38922, 39072, 63297,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 88707, 0, 3,
                                                                       86691, 61995, 86859,
                                                                       39072, 39222, 63507,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 88987, 0, 3,
                                                                       86859, 62121, 87027,
                                                                       39222, 39372, 63717,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 89267, 0, 3,
                                                                       87027, 62247, 87195,
                                                                       39372, 39522, 63927,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 89547, 0, 3,
                                                                       87195, 62373, 87363,
                                                                       39522, 39672, 64137,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 89827, 0, 3,
                                                                       87363, 62499, 87531,
                                                                       39672, 39822, 64347,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 90107, 0, 3,
                                                                       87531, 62625, 87699,
                                                                       39822, 39972, 64557,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 90387, 0, 3,
                                                                       87867, 62877, 88147,
                                                                       40272, 40497, 64767,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 90807, 0, 3,
                                                                       88147, 63087, 88427,
                                                                       40497, 40722, 65082,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 91227, 0, 3,
                                                                       88427, 63297, 88707,
                                                                       40722, 40947, 65397,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 91647, 0, 3,
                                                                       88707, 63507, 88987,
                                                                       40947, 41172, 65712,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 92067, 0, 3,
                                                                       88987, 63717, 89267,
                                                                       41172, 41397, 66027,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 92487, 0, 3,
                                                                       89267, 63927, 89547,
                                                                       41397, 41622, 66342,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 92907, 0, 3,
                                                                       89547, 64137, 89827,
                                                                       41622, 41847, 66657,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 93327, 0, 3,
                                                                       89827, 64347, 90107,
                                                                       41847, 42072, 66972,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 93747, 0, 3,
                                                                       90387, 64767, 90807,
                                                                       42522, 42837, 67287,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 94335, 0, 3,
                                                                       90807, 65082, 91227,
                                                                       42837, 43152, 67728,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 94923, 0, 3,
                                                                       91227, 65397, 91647,
                                                                       43152, 43467, 68169,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 95511, 0, 3,
                                                                       91647, 65712, 92067,
                                                                       43467, 43782, 68610,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 96099, 0, 3,
                                                                       92067, 66027, 92487,
                                                                       43782, 44097, 69051,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 96687, 0, 3,
                                                                       92487, 66342, 92907,
                                                                       44097, 44412, 69492,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 97275, 0, 3,
                                                                       92907, 66657, 93327,
                                                                       44412, 44727, 69933,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 97863, 0, 3,
                                                                       93747, 67287, 94335,
                                                                       45357, 45777, 70374,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 98647, 0, 3,
                                                                       94335, 67728, 94923,
                                                                       45777, 46197, 70962,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 99431, 0, 3,
                                                                       94923, 68169, 95511,
                                                                       46197, 46617, 71550,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 100215, 0, 3,
                                                                       95511, 68610, 96099,
                                                                       46617, 47037, 72138,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 100999, 0, 3,
                                                                       96099, 69051, 96687,
                                                                       47037, 47457, 72726,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 101783, 0, 3,
                                                                       96687, 69492, 97275,
                                                                       47457, 47877, 73314,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 102567, 0, 3,
                                                                       97863, 70374, 98647,
                                                                       48717, 49257, 73902,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 103575, 0, 3,
                                                                       98647, 70962, 99431,
                                                                       49257, 49797, 74658,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 104583, 0, 3,
                                                                       99431, 71550, 100215,
                                                                       49797, 50337, 75414,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 105591, 0, 3,
                                                                       100215, 72138, 100999,
                                                                       50337, 50877, 76170,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 106599, 0, 3,
                                                                       100999, 72726, 101783,
                                                                       50877, 51417, 76926,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 107607, 0, 3,
                                                                       102567, 73902, 103575,
                                                                       52497, 53172, 77682,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 108867, 0, 3,
                                                                       103575, 74658, 104583,
                                                                       53172, 53847, 78627,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 110127, 0, 3,
                                                                       104583, 75414, 105591,
                                                                       53847, 54522, 79572,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 111387, 0, 3,
                                                                       105591, 76170, 106599,
                                                                       54522, 55197, 80517,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 112647, 0, 3,
                                                                       107607, 77682, 108867,
                                                                       56547, 57372, 81462,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 114187, 0, 3,
                                                                       108867, 78627, 110127,
                                                                       57372, 58197, 82617,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 115727, 0, 3,
                                                                       110127, 79572, 111387,
                                                                       58197, 59022, 83772,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 117267, 3, 60672,
                                                                       60693, 84983, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 117303, 3, 60693,
                                                                       60714, 85011, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 117339, 3, 60714,
                                                                       60735, 85039, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 117375, 3, 60735,
                                                                       60756, 85067, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 117411, 3, 60756,
                                                                       60777, 85095, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 117447, 3, 60777,
                                                                       60798, 85123, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 117483, 3, 60798,
                                                                       60819, 85151, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 117519, 3, 60819,
                                                                       60840, 85179, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 117555, 3, 60840,
                                                                       60861, 85207, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 117591, 3, 60861,
                                                                       60882, 85235, ncols,
                                                                       gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 117627, 0, 3,
                                                                       117267, 84983, 117303,
                                                                       60924, 60987, 85431,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 117735, 0, 3,
                                                                       117303, 85011, 117339,
                                                                       60987, 61050, 85515,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 117843, 0, 3,
                                                                       117339, 85039, 117375,
                                                                       61050, 61113, 85599,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 117951, 0, 3,
                                                                       117375, 85067, 117411,
                                                                       61113, 61176, 85683,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 118059, 0, 3,
                                                                       117411, 85095, 117447,
                                                                       61176, 61239, 85767,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 118167, 0, 3,
                                                                       117447, 85123, 117483,
                                                                       61239, 61302, 85851,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 118275, 0, 3,
                                                                       117483, 85151, 117519,
                                                                       61302, 61365, 85935,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 118383, 0, 3,
                                                                       117519, 85179, 117555,
                                                                       61365, 61428, 86019,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 118491, 0, 3,
                                                                       117555, 85207, 117591,
                                                                       61428, 61491, 86103,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 118599, 0, 3,
                                                                       117627, 85431, 117735,
                                                                       61617, 61743, 86523,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 118815, 0, 3,
                                                                       117735, 85515, 117843,
                                                                       61743, 61869, 86691,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 119031, 0, 3,
                                                                       117843, 85599, 117951,
                                                                       61869, 61995, 86859,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 119247, 0, 3,
                                                                       117951, 85683, 118059,
                                                                       61995, 62121, 87027,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 119463, 0, 3,
                                                                       118059, 85767, 118167,
                                                                       62121, 62247, 87195,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 119679, 0, 3,
                                                                       118167, 85851, 118275,
                                                                       62247, 62373, 87363,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 119895, 0, 3,
                                                                       118275, 85935, 118383,
                                                                       62373, 62499, 87531,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 120111, 0, 3,
                                                                       118383, 86019, 118491,
                                                                       62499, 62625, 87699,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 120327, 0, 3,
                                                                       118599, 86523, 118815,
                                                                       62877, 63087, 88427,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 120687, 0, 3,
                                                                       118815, 86691, 119031,
                                                                       63087, 63297, 88707,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 121047, 0, 3,
                                                                       119031, 86859, 119247,
                                                                       63297, 63507, 88987,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 121407, 0, 3,
                                                                       119247, 87027, 119463,
                                                                       63507, 63717, 89267,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 121767, 0, 3,
                                                                       119463, 87195, 119679,
                                                                       63717, 63927, 89547,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 122127, 0, 3,
                                                                       119679, 87363, 119895,
                                                                       63927, 64137, 89827,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 122487, 0, 3,
                                                                       119895, 87531, 120111,
                                                                       64137, 64347, 90107,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 122847, 0, 3,
                                                                       120327, 88427, 120687,
                                                                       64767, 65082, 91227,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 123387, 0, 3,
                                                                       120687, 88707, 121047,
                                                                       65082, 65397, 91647,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 123927, 0, 3,
                                                                       121047, 88987, 121407,
                                                                       65397, 65712, 92067,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 124467, 0, 3,
                                                                       121407, 89267, 121767,
                                                                       65712, 66027, 92487,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 125007, 0, 3,
                                                                       121767, 89547, 122127,
                                                                       66027, 66342, 92907,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 125547, 0, 3,
                                                                       122127, 89827, 122487,
                                                                       66342, 66657, 93327,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 126087, 0, 3,
                                                                       122847, 91227, 123387,
                                                                       67287, 67728, 94923,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 126843, 0, 3,
                                                                       123387, 91647, 123927,
                                                                       67728, 68169, 95511,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 127599, 0, 3,
                                                                       123927, 92067, 124467,
                                                                       68169, 68610, 96099,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 128355, 0, 3,
                                                                       124467, 92487, 125007,
                                                                       68610, 69051, 96687,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 129111, 0, 3,
                                                                       125007, 92907, 125547,
                                                                       69051, 69492, 97275,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 129867, 0, 3,
                                                                       126087, 94923, 126843,
                                                                       70374, 70962, 99431,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 130875, 0, 3,
                                                                       126843, 95511, 127599,
                                                                       70962, 71550, 100215,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 131883, 0, 3,
                                                                       127599, 96099, 128355,
                                                                       71550, 72138, 100999,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 132891, 0, 3,
                                                                       128355, 96687, 129111,
                                                                       72138, 72726, 101783,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 133899, 0, 3,
                                                                       129867, 99431, 130875,
                                                                       73902, 74658, 104583,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 135195, 0, 3,
                                                                       130875, 100215, 131883,
                                                                       74658, 75414, 105591,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 136491, 0, 3,
                                                                       131883, 100999, 132891,
                                                                       75414, 76170, 106599,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 137787, 0, 3,
                                                                       133899, 104583, 135195,
                                                                       77682, 78627, 110127,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 139407, 0, 3,
                                                                       135195, 105591, 136491,
                                                                       78627, 79572, 111387,
                                                                       ncols, gamma, p, q);

                    compute_prim_smk_three_center_electron_repulsion_0(buffer, 141027, 0, 3,
                                                                       137787, 110127, 139407,
                                                                       81462, 82617, 115727,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 143007, 3, 84927,
                                                                       84955, 117267, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 143052, 3, 84955,
                                                                       84983, 117303, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 143097, 3, 84983,
                                                                       85011, 117339, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 143142, 3, 85011,
                                                                       85039, 117375, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 143187, 3, 85039,
                                                                       85067, 117411, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 143232, 3, 85067,
                                                                       85095, 117447, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 143277, 3, 85095,
                                                                       85123, 117483, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 143322, 3, 85123,
                                                                       85151, 117519, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 143367, 3, 85151,
                                                                       85179, 117555, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 143412, 3, 85179,
                                                                       85207, 117591, ncols,
                                                                       gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 143457, 0, 3,
                                                                       143007, 117267, 143052,
                                                                       85263, 85347, 117627,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 143592, 0, 3,
                                                                       143052, 117303, 143097,
                                                                       85347, 85431, 117735,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 143727, 0, 3,
                                                                       143097, 117339, 143142,
                                                                       85431, 85515, 117843,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 143862, 0, 3,
                                                                       143142, 117375, 143187,
                                                                       85515, 85599, 117951,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 143997, 0, 3,
                                                                       143187, 117411, 143232,
                                                                       85599, 85683, 118059,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 144132, 0, 3,
                                                                       143232, 117447, 143277,
                                                                       85683, 85767, 118167,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 144267, 0, 3,
                                                                       143277, 117483, 143322,
                                                                       85767, 85851, 118275,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 144402, 0, 3,
                                                                       143322, 117519, 143367,
                                                                       85851, 85935, 118383,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 144537, 0, 3,
                                                                       143367, 117555, 143412,
                                                                       85935, 86019, 118491,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 144672, 0, 3,
                                                                       143457, 117627, 143592,
                                                                       86187, 86355, 118599,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 144942, 0, 3,
                                                                       143592, 117735, 143727,
                                                                       86355, 86523, 118815,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 145212, 0, 3,
                                                                       143727, 117843, 143862,
                                                                       86523, 86691, 119031,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 145482, 0, 3,
                                                                       143862, 117951, 143997,
                                                                       86691, 86859, 119247,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 145752, 0, 3,
                                                                       143997, 118059, 144132,
                                                                       86859, 87027, 119463,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 146022, 0, 3,
                                                                       144132, 118167, 144267,
                                                                       87027, 87195, 119679,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 146292, 0, 3,
                                                                       144267, 118275, 144402,
                                                                       87195, 87363, 119895,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 146562, 0, 3,
                                                                       144402, 118383, 144537,
                                                                       87363, 87531, 120111,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 146832, 0, 3,
                                                                       144672, 118599, 144942,
                                                                       87867, 88147, 120327,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 147282, 0, 3,
                                                                       144942, 118815, 145212,
                                                                       88147, 88427, 120687,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 147732, 0, 3,
                                                                       145212, 119031, 145482,
                                                                       88427, 88707, 121047,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 148182, 0, 3,
                                                                       145482, 119247, 145752,
                                                                       88707, 88987, 121407,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 148632, 0, 3,
                                                                       145752, 119463, 146022,
                                                                       88987, 89267, 121767,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 149082, 0, 3,
                                                                       146022, 119679, 146292,
                                                                       89267, 89547, 122127,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 149532, 0, 3,
                                                                       146292, 119895, 146562,
                                                                       89547, 89827, 122487,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 149982, 0, 3,
                                                                       146832, 120327, 147282,
                                                                       90387, 90807, 122847,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 150657, 0, 3,
                                                                       147282, 120687, 147732,
                                                                       90807, 91227, 123387,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 151332, 0, 3,
                                                                       147732, 121047, 148182,
                                                                       91227, 91647, 123927,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 152007, 0, 3,
                                                                       148182, 121407, 148632,
                                                                       91647, 92067, 124467,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 152682, 0, 3,
                                                                       148632, 121767, 149082,
                                                                       92067, 92487, 125007,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 153357, 0, 3,
                                                                       149082, 122127, 149532,
                                                                       92487, 92907, 125547,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 154032, 0, 3,
                                                                       149982, 122847, 150657,
                                                                       93747, 94335, 126087,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 154977, 0, 3,
                                                                       150657, 123387, 151332,
                                                                       94335, 94923, 126843,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 155922, 0, 3,
                                                                       151332, 123927, 152007,
                                                                       94923, 95511, 127599,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 156867, 0, 3,
                                                                       152007, 124467, 152682,
                                                                       95511, 96099, 128355,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 157812, 0, 3,
                                                                       152682, 125007, 153357,
                                                                       96099, 96687, 129111,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 158757, 0, 3,
                                                                       154032, 126087, 154977,
                                                                       97863, 98647, 129867,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 160017, 0, 3,
                                                                       154977, 126843, 155922,
                                                                       98647, 99431, 130875,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 161277, 0, 3,
                                                                       155922, 127599, 156867,
                                                                       99431, 100215, 131883,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 162537, 0, 3,
                                                                       156867, 128355, 157812,
                                                                       100215, 100999, 132891,
                                                                       ncols, gamma, p, q);

                    compute_prim_skl_three_center_electron_repulsion_0(buffer, 163797, 0, 3,
                                                                       158757, 129867, 160017,
                                                                       102567, 103575, 133899,
                                                                       ncols, gamma, p, q);

                    compute_prim_skl_three_center_electron_repulsion_0(buffer, 165417, 0, 3,
                                                                       160017, 130875, 161277,
                                                                       103575, 104583, 135195,
                                                                       ncols, gamma, p, q);

                    compute_prim_skl_three_center_electron_repulsion_0(buffer, 167037, 0, 3,
                                                                       161277, 131883, 162537,
                                                                       104583, 105591, 136491,
                                                                       ncols, gamma, p, q);

                    compute_prim_sll_three_center_electron_repulsion_0(buffer, 168657, 0, 3,
                                                                       163797, 133899, 165417,
                                                                       107607, 108867, 137787,
                                                                       ncols, gamma, p, q);

                    compute_prim_sll_three_center_electron_repulsion_0(buffer, 170682, 0, 3,
                                                                       165417, 135195, 167037,
                                                                       108867, 110127, 139407,
                                                                       ncols, gamma, p, q);

                    compute_prim_sml_three_center_electron_repulsion_0(buffer, 172707, 0, 3,
                                                                       168657, 137787, 170682,
                                                                       112647, 114187, 141027,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 175182, 158757, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 176918, 163797, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 179150, 168657, 2025, ncols);

                    simdfunc::contract_primitives(buffer, 181940, 172707, 2475, ncols);
                }
            }
        }

        simdtrf::transform_l_inner(buffer, 176442, 175182, 28, 1, nmax);

        simdtrf::transform_l_inner(buffer, 178538, 176918, 36, 1, nmax);

        simdtrf::transform_l_inner(buffer, 181175, 179150, 45, 1, nmax);

        simdtrf::transform_l_inner(buffer, 184415, 181940, 55, 1, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 185350, 176442, 178538, 17, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 186778, 178538, 181175, 17, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 188614, 181175, 184415, 17, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 190909, 185350, 186778, 17, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 193765, 186778, 188614, 17, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 197437, 190909, 193765, 17, nmax);

        simdtrf::transform_i_inner(buffer, 202197, 197437, 10, 17, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 202197, 221, nmax);
    }

    for (size_t m = 0; m < 1547; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
