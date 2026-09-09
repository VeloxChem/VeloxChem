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


#include "SimdThreeCenterElectronRepulsionRecFHK.hpp"

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
compute_fhk_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_fhk_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 105454, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1155 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 105454, 85789, 5955, dimensions);

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
                                                        5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15},
                                                        ncols, fj, mu, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 22, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 25, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 28, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 31, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 34, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 37, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 40, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 43, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 46, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 49, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 52, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 55, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 58, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 61, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 64, 0, 3, 7, 8,
                                                                       22, 25, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 70, 0, 3, 8, 9,
                                                                       25, 28, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 76, 0, 3, 9, 10,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 82, 0, 3, 10, 11,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 88, 0, 3, 11, 12,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 94, 0, 3, 12, 13,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 100, 0, 3, 13, 14,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 106, 0, 3, 14, 15,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 112, 0, 3, 15, 16,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 118, 0, 3, 16, 17,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 124, 0, 3, 17, 18,
                                                                       52, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 130, 0, 3, 18, 19,
                                                                       55, 58, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 136, 0, 3, 19, 20,
                                                                       58, 61, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 142, 0, 3, 22, 25,
                                                                       64, 70, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 152, 0, 3, 25, 28,
                                                                       70, 76, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 162, 0, 3, 28, 31,
                                                                       76, 82, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 172, 0, 3, 31, 34,
                                                                       82, 88, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 182, 0, 3, 34, 37,
                                                                       88, 94, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 192, 0, 3, 37, 40,
                                                                       94, 100, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 202, 0, 3, 40, 43,
                                                                       100, 106, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 212, 0, 3, 43, 46,
                                                                       106, 112, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 222, 0, 3, 46, 49,
                                                                       112, 118, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 232, 0, 3, 49, 52,
                                                                       118, 124, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 242, 0, 3, 52, 55,
                                                                       124, 130, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 252, 0, 3, 55, 58,
                                                                       130, 136, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 262, 0, 3, 64, 70,
                                                                       142, 152, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 277, 0, 3, 70, 76,
                                                                       152, 162, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 292, 0, 3, 76, 82,
                                                                       162, 172, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 307, 0, 3, 82, 88,
                                                                       172, 182, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 322, 0, 3, 88, 94,
                                                                       182, 192, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 337, 0, 3, 94,
                                                                       100, 192, 202, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 352, 0, 3, 100,
                                                                       106, 202, 212, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 367, 0, 3, 106,
                                                                       112, 212, 222, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 382, 0, 3, 112,
                                                                       118, 222, 232, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 397, 0, 3, 118,
                                                                       124, 232, 242, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 412, 0, 3, 124,
                                                                       130, 242, 252, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 427, 0, 3, 142,
                                                                       152, 262, 277, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 448, 0, 3, 152,
                                                                       162, 277, 292, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 469, 0, 3, 162,
                                                                       172, 292, 307, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 490, 0, 3, 172,
                                                                       182, 307, 322, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 511, 0, 3, 182,
                                                                       192, 322, 337, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 532, 0, 3, 192,
                                                                       202, 337, 352, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 553, 0, 3, 202,
                                                                       212, 352, 367, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 574, 0, 3, 212,
                                                                       222, 367, 382, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 595, 0, 3, 222,
                                                                       232, 382, 397, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 616, 0, 3, 232,
                                                                       242, 397, 412, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 637, 0, 3, 262,
                                                                       277, 427, 448, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 665, 0, 3, 277,
                                                                       292, 448, 469, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 693, 0, 3, 292,
                                                                       307, 469, 490, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 721, 0, 3, 307,
                                                                       322, 490, 511, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 749, 0, 3, 322,
                                                                       337, 511, 532, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 777, 0, 3, 337,
                                                                       352, 532, 553, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 805, 0, 3, 352,
                                                                       367, 553, 574, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 833, 0, 3, 367,
                                                                       382, 574, 595, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 861, 0, 3, 382,
                                                                       397, 595, 616, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 889, 0, 3, 427,
                                                                       448, 637, 665, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 925, 0, 3, 448,
                                                                       469, 665, 693, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 961, 0, 3, 469,
                                                                       490, 693, 721, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 997, 0, 3, 490,
                                                                       511, 721, 749, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1033, 0, 3, 511,
                                                                       532, 749, 777, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1069, 0, 3, 532,
                                                                       553, 777, 805, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1105, 0, 3, 553,
                                                                       574, 805, 833, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1141, 0, 3, 574,
                                                                       595, 833, 861, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1177, 0, 3, 637,
                                                                       665, 889, 925, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1222, 0, 3, 665,
                                                                       693, 925, 961, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1267, 0, 3, 693,
                                                                       721, 961, 997, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1312, 0, 3, 721,
                                                                       749, 997, 1033, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1357, 0, 3, 749,
                                                                       777, 1033, 1069, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1402, 0, 3, 777,
                                                                       805, 1069, 1105, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1447, 0, 3, 805,
                                                                       833, 1105, 1141, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1492, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1495, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1498, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1501, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1504, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1507, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1510, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1513, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1516, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1519, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1522, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1525, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1528, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1531, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1534, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1537, 3, 9, 28,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1546, 3, 10, 31,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1555, 3, 11, 34,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1564, 3, 12, 37,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1573, 3, 13, 40,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1582, 3, 14, 43,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1591, 3, 15, 46,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1600, 3, 16, 49,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1609, 3, 17, 52,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1618, 3, 18, 55,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1627, 3, 19, 58,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1636, 3, 20, 61,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1645, 3, 22, 64,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1663, 3, 25, 70,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1681, 3, 28, 76,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1699, 3, 31, 82,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1717, 3, 34, 88,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1735, 3, 37, 94,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1753, 3, 40, 100,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1771, 3, 43, 106,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1789, 3, 46, 112,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1807, 3, 49, 118,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1825, 3, 52, 124,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1843, 3, 55, 130,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1861, 3, 58, 136,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1879, 3, 64, 142,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1909, 3, 70, 152,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1939, 3, 76, 162,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1969, 3, 82, 172,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1999, 3, 88, 182,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2029, 3, 94, 192,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2059, 3, 100, 202,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2089, 3, 106, 212,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2119, 3, 112, 222,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2149, 3, 118, 232,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2179, 3, 124, 242,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2209, 3, 130, 252,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2239, 3, 142, 262,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2284, 3, 152, 277,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2329, 3, 162, 292,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2374, 3, 172, 307,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2419, 3, 182, 322,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2464, 3, 192, 337,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2509, 3, 202, 352,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2554, 3, 212, 367,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2599, 3, 222, 382,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2644, 3, 232, 397,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2689, 3, 242, 412,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2734, 3, 262, 427,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2797, 3, 277, 448,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2860, 3, 292, 469,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2923, 3, 307, 490,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2986, 3, 322, 511,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3049, 3, 337, 532,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3112, 3, 352, 553,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3175, 3, 367, 574,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3238, 3, 382, 595,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3301, 3, 397, 616,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3364, 3, 427, 637,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3448, 3, 448, 665,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3532, 3, 469, 693,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3616, 3, 490, 721,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3700, 3, 511, 749,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3784, 3, 532, 777,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3868, 3, 553, 805,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3952, 3, 574, 833,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4036, 3, 595, 861,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4120, 3, 637, 889,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4228, 3, 665, 925,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4336, 3, 693, 961,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4444, 3, 721, 997,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4552, 3, 749,
                                                                       1033, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4660, 3, 777,
                                                                       1069, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4768, 3, 805,
                                                                       1105, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4876, 3, 833,
                                                                       1141, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 4984, 3, 889,
                                                                       1177, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5119, 3, 925,
                                                                       1222, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5254, 3, 961,
                                                                       1267, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5389, 3, 997,
                                                                       1312, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5524, 3, 1033,
                                                                       1357, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5659, 3, 1069,
                                                                       1402, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5794, 3, 1105,
                                                                       1447, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5929, 3, 7, 8,
                                                                       1498, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5935, 3, 8, 9,
                                                                       1501, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5941, 3, 9, 10,
                                                                       1504, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5947, 3, 10, 11,
                                                                       1507, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5953, 3, 11, 12,
                                                                       1510, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5959, 3, 12, 13,
                                                                       1513, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5965, 3, 13, 14,
                                                                       1516, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5971, 3, 14, 15,
                                                                       1519, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5977, 3, 15, 16,
                                                                       1522, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5983, 3, 16, 17,
                                                                       1525, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5989, 3, 17, 18,
                                                                       1528, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5995, 3, 18, 19,
                                                                       1531, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6001, 3, 19, 20,
                                                                       1534, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6007, 0, 3, 5929,
                                                                       1498, 5935, 1537, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6025, 0, 3, 5935,
                                                                       1501, 5941, 1546, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6043, 0, 3, 5941,
                                                                       1504, 5947, 1555, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6061, 0, 3, 5947,
                                                                       1507, 5953, 1564, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6079, 0, 3, 5953,
                                                                       1510, 5959, 1573, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6097, 0, 3, 5959,
                                                                       1513, 5965, 1582, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6115, 0, 3, 5965,
                                                                       1516, 5971, 1591, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6133, 0, 3, 5971,
                                                                       1519, 5977, 1600, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6151, 0, 3, 5977,
                                                                       1522, 5983, 1609, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6169, 0, 3, 5983,
                                                                       1525, 5989, 1618, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6187, 0, 3, 5989,
                                                                       1528, 5995, 1627, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6205, 0, 3, 5995,
                                                                       1531, 6001, 1636, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6223, 0, 3, 6007,
                                                                       1537, 6025, 64, 70, 1681,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6259, 0, 3, 6025,
                                                                       1546, 6043, 70, 76, 1699,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6295, 0, 3, 6043,
                                                                       1555, 6061, 76, 82, 1717,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6331, 0, 3, 6061,
                                                                       1564, 6079, 82, 88, 1735,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6367, 0, 3, 6079,
                                                                       1573, 6097, 88, 94, 1753,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6403, 0, 3, 6097,
                                                                       1582, 6115, 94, 100, 1771,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6439, 0, 3, 6115,
                                                                       1591, 6133, 100, 106,
                                                                       1789, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6475, 0, 3, 6133,
                                                                       1600, 6151, 106, 112,
                                                                       1807, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6511, 0, 3, 6151,
                                                                       1609, 6169, 112, 118,
                                                                       1825, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6547, 0, 3, 6169,
                                                                       1618, 6187, 118, 124,
                                                                       1843, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6583, 0, 3, 6187,
                                                                       1627, 6205, 124, 130,
                                                                       1861, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6619, 0, 3, 6223,
                                                                       1681, 6259, 142, 152,
                                                                       1939, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6679, 0, 3, 6259,
                                                                       1699, 6295, 152, 162,
                                                                       1969, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6739, 0, 3, 6295,
                                                                       1717, 6331, 162, 172,
                                                                       1999, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6799, 0, 3, 6331,
                                                                       1735, 6367, 172, 182,
                                                                       2029, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6859, 0, 3, 6367,
                                                                       1753, 6403, 182, 192,
                                                                       2059, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6919, 0, 3, 6403,
                                                                       1771, 6439, 192, 202,
                                                                       2089, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6979, 0, 3, 6439,
                                                                       1789, 6475, 202, 212,
                                                                       2119, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7039, 0, 3, 6475,
                                                                       1807, 6511, 212, 222,
                                                                       2149, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7099, 0, 3, 6511,
                                                                       1825, 6547, 222, 232,
                                                                       2179, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7159, 0, 3, 6547,
                                                                       1843, 6583, 232, 242,
                                                                       2209, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7219, 0, 3, 6619,
                                                                       1939, 6679, 262, 277,
                                                                       2329, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7309, 0, 3, 6679,
                                                                       1969, 6739, 277, 292,
                                                                       2374, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7399, 0, 3, 6739,
                                                                       1999, 6799, 292, 307,
                                                                       2419, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7489, 0, 3, 6799,
                                                                       2029, 6859, 307, 322,
                                                                       2464, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7579, 0, 3, 6859,
                                                                       2059, 6919, 322, 337,
                                                                       2509, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7669, 0, 3, 6919,
                                                                       2089, 6979, 337, 352,
                                                                       2554, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7759, 0, 3, 6979,
                                                                       2119, 7039, 352, 367,
                                                                       2599, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7849, 0, 3, 7039,
                                                                       2149, 7099, 367, 382,
                                                                       2644, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7939, 0, 3, 7099,
                                                                       2179, 7159, 382, 397,
                                                                       2689, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8029, 0, 3, 7219,
                                                                       2329, 7309, 427, 448,
                                                                       2860, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8155, 0, 3, 7309,
                                                                       2374, 7399, 448, 469,
                                                                       2923, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8281, 0, 3, 7399,
                                                                       2419, 7489, 469, 490,
                                                                       2986, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8407, 0, 3, 7489,
                                                                       2464, 7579, 490, 511,
                                                                       3049, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8533, 0, 3, 7579,
                                                                       2509, 7669, 511, 532,
                                                                       3112, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8659, 0, 3, 7669,
                                                                       2554, 7759, 532, 553,
                                                                       3175, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8785, 0, 3, 7759,
                                                                       2599, 7849, 553, 574,
                                                                       3238, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8911, 0, 3, 7849,
                                                                       2644, 7939, 574, 595,
                                                                       3301, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 9037, 0, 3, 8029,
                                                                       2860, 8155, 637, 665,
                                                                       3532, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 9205, 0, 3, 8155,
                                                                       2923, 8281, 665, 693,
                                                                       3616, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 9373, 0, 3, 8281,
                                                                       2986, 8407, 693, 721,
                                                                       3700, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 9541, 0, 3, 8407,
                                                                       3049, 8533, 721, 749,
                                                                       3784, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 9709, 0, 3, 8533,
                                                                       3112, 8659, 749, 777,
                                                                       3868, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 9877, 0, 3, 8659,
                                                                       3175, 8785, 777, 805,
                                                                       3952, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 10045, 0, 3, 8785,
                                                                       3238, 8911, 805, 833,
                                                                       4036, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 10213, 0, 3, 9037,
                                                                       3532, 9205, 889, 925,
                                                                       4336, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 10429, 0, 3, 9205,
                                                                       3616, 9373, 925, 961,
                                                                       4444, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 10645, 0, 3, 9373,
                                                                       3700, 9541, 961, 997,
                                                                       4552, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 10861, 0, 3, 9541,
                                                                       3784, 9709, 997, 1033,
                                                                       4660, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 11077, 0, 3, 9709,
                                                                       3868, 9877, 1033, 1069,
                                                                       4768, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 11293, 0, 3, 9877,
                                                                       3952, 10045, 1069, 1105,
                                                                       4876, ncols, gamma, p,
                                                                       q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 11509, 0, 3,
                                                                       10213, 4336, 10429, 1177,
                                                                       1222, 5254, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 11779, 0, 3,
                                                                       10429, 4444, 10645, 1222,
                                                                       1267, 5389, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 12049, 0, 3,
                                                                       10645, 4552, 10861, 1267,
                                                                       1312, 5524, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 12319, 0, 3,
                                                                       10861, 4660, 11077, 1312,
                                                                       1357, 5659, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 12589, 0, 3,
                                                                       11077, 4768, 11293, 1357,
                                                                       1402, 5794, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12859, 3, 1492,
                                                                       1495, 5929, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12869, 3, 1495,
                                                                       1498, 5935, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12879, 3, 1498,
                                                                       1501, 5941, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12889, 3, 1501,
                                                                       1504, 5947, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12899, 3, 1504,
                                                                       1507, 5953, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12909, 3, 1507,
                                                                       1510, 5959, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12919, 3, 1510,
                                                                       1513, 5965, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12929, 3, 1513,
                                                                       1516, 5971, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12939, 3, 1516,
                                                                       1519, 5977, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12949, 3, 1519,
                                                                       1522, 5983, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12959, 3, 1522,
                                                                       1525, 5989, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12969, 3, 1525,
                                                                       1528, 5995, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12979, 3, 1528,
                                                                       1531, 6001, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12989, 0, 3,
                                                                       12859, 5929, 12869, 6007,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13019, 0, 3,
                                                                       12869, 5935, 12879, 6025,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13049, 0, 3,
                                                                       12879, 5941, 12889, 6043,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13079, 0, 3,
                                                                       12889, 5947, 12899, 6061,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13109, 0, 3,
                                                                       12899, 5953, 12909, 6079,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13139, 0, 3,
                                                                       12909, 5959, 12919, 6097,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13169, 0, 3,
                                                                       12919, 5965, 12929, 6115,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13199, 0, 3,
                                                                       12929, 5971, 12939, 6133,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13229, 0, 3,
                                                                       12939, 5977, 12949, 6151,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13259, 0, 3,
                                                                       12949, 5983, 12959, 6169,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13289, 0, 3,
                                                                       12959, 5989, 12969, 6187,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13319, 0, 3,
                                                                       12969, 5995, 12979, 6205,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13349, 0, 3,
                                                                       12989, 6007, 13019, 1645,
                                                                       1663, 6223, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13409, 0, 3,
                                                                       13019, 6025, 13049, 1663,
                                                                       1681, 6259, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13469, 0, 3,
                                                                       13049, 6043, 13079, 1681,
                                                                       1699, 6295, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13529, 0, 3,
                                                                       13079, 6061, 13109, 1699,
                                                                       1717, 6331, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13589, 0, 3,
                                                                       13109, 6079, 13139, 1717,
                                                                       1735, 6367, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13649, 0, 3,
                                                                       13139, 6097, 13169, 1735,
                                                                       1753, 6403, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13709, 0, 3,
                                                                       13169, 6115, 13199, 1753,
                                                                       1771, 6439, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13769, 0, 3,
                                                                       13199, 6133, 13229, 1771,
                                                                       1789, 6475, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13829, 0, 3,
                                                                       13229, 6151, 13259, 1789,
                                                                       1807, 6511, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13889, 0, 3,
                                                                       13259, 6169, 13289, 1807,
                                                                       1825, 6547, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13949, 0, 3,
                                                                       13289, 6187, 13319, 1825,
                                                                       1843, 6583, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14009, 0, 3,
                                                                       13349, 6223, 13409, 1879,
                                                                       1909, 6619, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14109, 0, 3,
                                                                       13409, 6259, 13469, 1909,
                                                                       1939, 6679, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14209, 0, 3,
                                                                       13469, 6295, 13529, 1939,
                                                                       1969, 6739, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14309, 0, 3,
                                                                       13529, 6331, 13589, 1969,
                                                                       1999, 6799, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14409, 0, 3,
                                                                       13589, 6367, 13649, 1999,
                                                                       2029, 6859, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14509, 0, 3,
                                                                       13649, 6403, 13709, 2029,
                                                                       2059, 6919, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14609, 0, 3,
                                                                       13709, 6439, 13769, 2059,
                                                                       2089, 6979, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14709, 0, 3,
                                                                       13769, 6475, 13829, 2089,
                                                                       2119, 7039, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14809, 0, 3,
                                                                       13829, 6511, 13889, 2119,
                                                                       2149, 7099, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14909, 0, 3,
                                                                       13889, 6547, 13949, 2149,
                                                                       2179, 7159, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 15009, 0, 3,
                                                                       14009, 6619, 14109, 2239,
                                                                       2284, 7219, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 15159, 0, 3,
                                                                       14109, 6679, 14209, 2284,
                                                                       2329, 7309, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 15309, 0, 3,
                                                                       14209, 6739, 14309, 2329,
                                                                       2374, 7399, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 15459, 0, 3,
                                                                       14309, 6799, 14409, 2374,
                                                                       2419, 7489, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 15609, 0, 3,
                                                                       14409, 6859, 14509, 2419,
                                                                       2464, 7579, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 15759, 0, 3,
                                                                       14509, 6919, 14609, 2464,
                                                                       2509, 7669, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 15909, 0, 3,
                                                                       14609, 6979, 14709, 2509,
                                                                       2554, 7759, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 16059, 0, 3,
                                                                       14709, 7039, 14809, 2554,
                                                                       2599, 7849, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 16209, 0, 3,
                                                                       14809, 7099, 14909, 2599,
                                                                       2644, 7939, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 16359, 0, 3,
                                                                       15009, 7219, 15159, 2734,
                                                                       2797, 8029, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 16569, 0, 3,
                                                                       15159, 7309, 15309, 2797,
                                                                       2860, 8155, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 16779, 0, 3,
                                                                       15309, 7399, 15459, 2860,
                                                                       2923, 8281, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 16989, 0, 3,
                                                                       15459, 7489, 15609, 2923,
                                                                       2986, 8407, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 17199, 0, 3,
                                                                       15609, 7579, 15759, 2986,
                                                                       3049, 8533, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 17409, 0, 3,
                                                                       15759, 7669, 15909, 3049,
                                                                       3112, 8659, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 17619, 0, 3,
                                                                       15909, 7759, 16059, 3112,
                                                                       3175, 8785, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 17829, 0, 3,
                                                                       16059, 7849, 16209, 3175,
                                                                       3238, 8911, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 18039, 0, 3,
                                                                       16359, 8029, 16569, 3364,
                                                                       3448, 9037, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 18319, 0, 3,
                                                                       16569, 8155, 16779, 3448,
                                                                       3532, 9205, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 18599, 0, 3,
                                                                       16779, 8281, 16989, 3532,
                                                                       3616, 9373, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 18879, 0, 3,
                                                                       16989, 8407, 17199, 3616,
                                                                       3700, 9541, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 19159, 0, 3,
                                                                       17199, 8533, 17409, 3700,
                                                                       3784, 9709, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 19439, 0, 3,
                                                                       17409, 8659, 17619, 3784,
                                                                       3868, 9877, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 19719, 0, 3,
                                                                       17619, 8785, 17829, 3868,
                                                                       3952, 10045, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 19999, 0, 3,
                                                                       18039, 9037, 18319, 4120,
                                                                       4228, 10213, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 20359, 0, 3,
                                                                       18319, 9205, 18599, 4228,
                                                                       4336, 10429, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 20719, 0, 3,
                                                                       18599, 9373, 18879, 4336,
                                                                       4444, 10645, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 21079, 0, 3,
                                                                       18879, 9541, 19159, 4444,
                                                                       4552, 10861, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 21439, 0, 3,
                                                                       19159, 9709, 19439, 4552,
                                                                       4660, 11077, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 21799, 0, 3,
                                                                       19439, 9877, 19719, 4660,
                                                                       4768, 11293, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 22159, 0, 3,
                                                                       19999, 10213, 20359, 4984,
                                                                       5119, 11509, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 22609, 0, 3,
                                                                       20359, 10429, 20719, 5119,
                                                                       5254, 11779, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 23059, 0, 3,
                                                                       20719, 10645, 21079, 5254,
                                                                       5389, 12049, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 23509, 0, 3,
                                                                       21079, 10861, 21439, 5389,
                                                                       5524, 12319, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 23959, 0, 3,
                                                                       21439, 11077, 21799, 5524,
                                                                       5659, 12589, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24409, 3, 5929,
                                                                       5935, 12879, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24424, 3, 5935,
                                                                       5941, 12889, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24439, 3, 5941,
                                                                       5947, 12899, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24454, 3, 5947,
                                                                       5953, 12909, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24469, 3, 5953,
                                                                       5959, 12919, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24484, 3, 5959,
                                                                       5965, 12929, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24499, 3, 5965,
                                                                       5971, 12939, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24514, 3, 5971,
                                                                       5977, 12949, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24529, 3, 5977,
                                                                       5983, 12959, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24544, 3, 5983,
                                                                       5989, 12969, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24559, 3, 5989,
                                                                       5995, 12979, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 24574, 0, 3,
                                                                       24409, 12879, 24424, 6007,
                                                                       6025, 13049, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 24619, 0, 3,
                                                                       24424, 12889, 24439, 6025,
                                                                       6043, 13079, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 24664, 0, 3,
                                                                       24439, 12899, 24454, 6043,
                                                                       6061, 13109, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 24709, 0, 3,
                                                                       24454, 12909, 24469, 6061,
                                                                       6079, 13139, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 24754, 0, 3,
                                                                       24469, 12919, 24484, 6079,
                                                                       6097, 13169, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 24799, 0, 3,
                                                                       24484, 12929, 24499, 6097,
                                                                       6115, 13199, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 24844, 0, 3,
                                                                       24499, 12939, 24514, 6115,
                                                                       6133, 13229, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 24889, 0, 3,
                                                                       24514, 12949, 24529, 6133,
                                                                       6151, 13259, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 24934, 0, 3,
                                                                       24529, 12959, 24544, 6151,
                                                                       6169, 13289, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 24979, 0, 3,
                                                                       24544, 12969, 24559, 6169,
                                                                       6187, 13319, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 25024, 0, 3,
                                                                       24574, 13049, 24619, 6223,
                                                                       6259, 13469, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 25114, 0, 3,
                                                                       24619, 13079, 24664, 6259,
                                                                       6295, 13529, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 25204, 0, 3,
                                                                       24664, 13109, 24709, 6295,
                                                                       6331, 13589, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 25294, 0, 3,
                                                                       24709, 13139, 24754, 6331,
                                                                       6367, 13649, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 25384, 0, 3,
                                                                       24754, 13169, 24799, 6367,
                                                                       6403, 13709, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 25474, 0, 3,
                                                                       24799, 13199, 24844, 6403,
                                                                       6439, 13769, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 25564, 0, 3,
                                                                       24844, 13229, 24889, 6439,
                                                                       6475, 13829, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 25654, 0, 3,
                                                                       24889, 13259, 24934, 6475,
                                                                       6511, 13889, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 25744, 0, 3,
                                                                       24934, 13289, 24979, 6511,
                                                                       6547, 13949, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 25834, 0, 3,
                                                                       25024, 13469, 25114, 6619,
                                                                       6679, 14209, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 25984, 0, 3,
                                                                       25114, 13529, 25204, 6679,
                                                                       6739, 14309, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 26134, 0, 3,
                                                                       25204, 13589, 25294, 6739,
                                                                       6799, 14409, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 26284, 0, 3,
                                                                       25294, 13649, 25384, 6799,
                                                                       6859, 14509, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 26434, 0, 3,
                                                                       25384, 13709, 25474, 6859,
                                                                       6919, 14609, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 26584, 0, 3,
                                                                       25474, 13769, 25564, 6919,
                                                                       6979, 14709, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 26734, 0, 3,
                                                                       25564, 13829, 25654, 6979,
                                                                       7039, 14809, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 26884, 0, 3,
                                                                       25654, 13889, 25744, 7039,
                                                                       7099, 14909, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 27034, 0, 3,
                                                                       25834, 14209, 25984, 7219,
                                                                       7309, 15309, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 27259, 0, 3,
                                                                       25984, 14309, 26134, 7309,
                                                                       7399, 15459, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 27484, 0, 3,
                                                                       26134, 14409, 26284, 7399,
                                                                       7489, 15609, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 27709, 0, 3,
                                                                       26284, 14509, 26434, 7489,
                                                                       7579, 15759, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 27934, 0, 3,
                                                                       26434, 14609, 26584, 7579,
                                                                       7669, 15909, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 28159, 0, 3,
                                                                       26584, 14709, 26734, 7669,
                                                                       7759, 16059, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 28384, 0, 3,
                                                                       26734, 14809, 26884, 7759,
                                                                       7849, 16209, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 28609, 0, 3,
                                                                       27034, 15309, 27259, 8029,
                                                                       8155, 16779, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 28924, 0, 3,
                                                                       27259, 15459, 27484, 8155,
                                                                       8281, 16989, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 29239, 0, 3,
                                                                       27484, 15609, 27709, 8281,
                                                                       8407, 17199, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 29554, 0, 3,
                                                                       27709, 15759, 27934, 8407,
                                                                       8533, 17409, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 29869, 0, 3,
                                                                       27934, 15909, 28159, 8533,
                                                                       8659, 17619, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 30184, 0, 3,
                                                                       28159, 16059, 28384, 8659,
                                                                       8785, 17829, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 30499, 0, 3,
                                                                       28609, 16779, 28924, 9037,
                                                                       9205, 18599, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 30919, 0, 3,
                                                                       28924, 16989, 29239, 9205,
                                                                       9373, 18879, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 31339, 0, 3,
                                                                       29239, 17199, 29554, 9373,
                                                                       9541, 19159, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 31759, 0, 3,
                                                                       29554, 17409, 29869, 9541,
                                                                       9709, 19439, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 32179, 0, 3,
                                                                       29869, 17619, 30184, 9709,
                                                                       9877, 19719, ncols, gamma,
                                                                       p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 32599, 0, 3,
                                                                       30499, 18599, 30919,
                                                                       10213, 10429, 20719,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 33139, 0, 3,
                                                                       30919, 18879, 31339,
                                                                       10429, 10645, 21079,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 33679, 0, 3,
                                                                       31339, 19159, 31759,
                                                                       10645, 10861, 21439,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 34219, 0, 3,
                                                                       31759, 19439, 32179,
                                                                       10861, 11077, 21799,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 34759, 0, 3,
                                                                       32599, 20719, 33139,
                                                                       11509, 11779, 23059,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 35434, 0, 3,
                                                                       33139, 21079, 33679,
                                                                       11779, 12049, 23509,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 36109, 0, 3,
                                                                       33679, 21439, 34219,
                                                                       12049, 12319, 23959,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 36784, 3, 12859,
                                                                       12869, 24409, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 36805, 3, 12869,
                                                                       12879, 24424, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 36826, 3, 12879,
                                                                       12889, 24439, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 36847, 3, 12889,
                                                                       12899, 24454, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 36868, 3, 12899,
                                                                       12909, 24469, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 36889, 3, 12909,
                                                                       12919, 24484, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 36910, 3, 12919,
                                                                       12929, 24499, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 36931, 3, 12929,
                                                                       12939, 24514, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 36952, 3, 12939,
                                                                       12949, 24529, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 36973, 3, 12949,
                                                                       12959, 24544, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 36994, 3, 12959,
                                                                       12969, 24559, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 37015, 0, 3,
                                                                       36784, 24409, 36805,
                                                                       12989, 13019, 24574,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 37078, 0, 3,
                                                                       36805, 24424, 36826,
                                                                       13019, 13049, 24619,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 37141, 0, 3,
                                                                       36826, 24439, 36847,
                                                                       13049, 13079, 24664,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 37204, 0, 3,
                                                                       36847, 24454, 36868,
                                                                       13079, 13109, 24709,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 37267, 0, 3,
                                                                       36868, 24469, 36889,
                                                                       13109, 13139, 24754,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 37330, 0, 3,
                                                                       36889, 24484, 36910,
                                                                       13139, 13169, 24799,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 37393, 0, 3,
                                                                       36910, 24499, 36931,
                                                                       13169, 13199, 24844,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 37456, 0, 3,
                                                                       36931, 24514, 36952,
                                                                       13199, 13229, 24889,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 37519, 0, 3,
                                                                       36952, 24529, 36973,
                                                                       13229, 13259, 24934,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 37582, 0, 3,
                                                                       36973, 24544, 36994,
                                                                       13259, 13289, 24979,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 37645, 0, 3,
                                                                       37015, 24574, 37078,
                                                                       13349, 13409, 25024,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 37771, 0, 3,
                                                                       37078, 24619, 37141,
                                                                       13409, 13469, 25114,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 37897, 0, 3,
                                                                       37141, 24664, 37204,
                                                                       13469, 13529, 25204,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 38023, 0, 3,
                                                                       37204, 24709, 37267,
                                                                       13529, 13589, 25294,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 38149, 0, 3,
                                                                       37267, 24754, 37330,
                                                                       13589, 13649, 25384,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 38275, 0, 3,
                                                                       37330, 24799, 37393,
                                                                       13649, 13709, 25474,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 38401, 0, 3,
                                                                       37393, 24844, 37456,
                                                                       13709, 13769, 25564,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 38527, 0, 3,
                                                                       37456, 24889, 37519,
                                                                       13769, 13829, 25654,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 38653, 0, 3,
                                                                       37519, 24934, 37582,
                                                                       13829, 13889, 25744,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 38779, 0, 3,
                                                                       37645, 25024, 37771,
                                                                       14009, 14109, 25834,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 38989, 0, 3,
                                                                       37771, 25114, 37897,
                                                                       14109, 14209, 25984,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 39199, 0, 3,
                                                                       37897, 25204, 38023,
                                                                       14209, 14309, 26134,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 39409, 0, 3,
                                                                       38023, 25294, 38149,
                                                                       14309, 14409, 26284,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 39619, 0, 3,
                                                                       38149, 25384, 38275,
                                                                       14409, 14509, 26434,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 39829, 0, 3,
                                                                       38275, 25474, 38401,
                                                                       14509, 14609, 26584,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 40039, 0, 3,
                                                                       38401, 25564, 38527,
                                                                       14609, 14709, 26734,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 40249, 0, 3,
                                                                       38527, 25654, 38653,
                                                                       14709, 14809, 26884,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 40459, 0, 3,
                                                                       38779, 25834, 38989,
                                                                       15009, 15159, 27034,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 40774, 0, 3,
                                                                       38989, 25984, 39199,
                                                                       15159, 15309, 27259,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 41089, 0, 3,
                                                                       39199, 26134, 39409,
                                                                       15309, 15459, 27484,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 41404, 0, 3,
                                                                       39409, 26284, 39619,
                                                                       15459, 15609, 27709,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 41719, 0, 3,
                                                                       39619, 26434, 39829,
                                                                       15609, 15759, 27934,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 42034, 0, 3,
                                                                       39829, 26584, 40039,
                                                                       15759, 15909, 28159,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 42349, 0, 3,
                                                                       40039, 26734, 40249,
                                                                       15909, 16059, 28384,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 42664, 0, 3,
                                                                       40459, 27034, 40774,
                                                                       16359, 16569, 28609,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 43105, 0, 3,
                                                                       40774, 27259, 41089,
                                                                       16569, 16779, 28924,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 43546, 0, 3,
                                                                       41089, 27484, 41404,
                                                                       16779, 16989, 29239,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 43987, 0, 3,
                                                                       41404, 27709, 41719,
                                                                       16989, 17199, 29554,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 44428, 0, 3,
                                                                       41719, 27934, 42034,
                                                                       17199, 17409, 29869,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 44869, 0, 3,
                                                                       42034, 28159, 42349,
                                                                       17409, 17619, 30184,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 45310, 0, 3,
                                                                       42664, 28609, 43105,
                                                                       18039, 18319, 30499,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 45898, 0, 3,
                                                                       43105, 28924, 43546,
                                                                       18319, 18599, 30919,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 46486, 0, 3,
                                                                       43546, 29239, 43987,
                                                                       18599, 18879, 31339,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 47074, 0, 3,
                                                                       43987, 29554, 44428,
                                                                       18879, 19159, 31759,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 47662, 0, 3,
                                                                       44428, 29869, 44869,
                                                                       19159, 19439, 32179,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 48250, 0, 3,
                                                                       45310, 30499, 45898,
                                                                       19999, 20359, 32599,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 49006, 0, 3,
                                                                       45898, 30919, 46486,
                                                                       20359, 20719, 33139,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 49762, 0, 3,
                                                                       46486, 31339, 47074,
                                                                       20719, 21079, 33679,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 50518, 0, 3,
                                                                       47074, 31759, 47662,
                                                                       21079, 21439, 34219,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 51274, 0, 3,
                                                                       48250, 32599, 49006,
                                                                       22159, 22609, 34759,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 52219, 0, 3,
                                                                       49006, 33139, 49762,
                                                                       22609, 23059, 35434,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 53164, 0, 3,
                                                                       49762, 33679, 50518,
                                                                       23059, 23509, 36109,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 54109, 3, 24409,
                                                                       24424, 36826, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 54137, 3, 24424,
                                                                       24439, 36847, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 54165, 3, 24439,
                                                                       24454, 36868, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 54193, 3, 24454,
                                                                       24469, 36889, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 54221, 3, 24469,
                                                                       24484, 36910, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 54249, 3, 24484,
                                                                       24499, 36931, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 54277, 3, 24499,
                                                                       24514, 36952, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 54305, 3, 24514,
                                                                       24529, 36973, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 54333, 3, 24529,
                                                                       24544, 36994, ncols,
                                                                       gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 54361, 0, 3,
                                                                       54109, 36826, 54137,
                                                                       24574, 24619, 37141,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 54445, 0, 3,
                                                                       54137, 36847, 54165,
                                                                       24619, 24664, 37204,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 54529, 0, 3,
                                                                       54165, 36868, 54193,
                                                                       24664, 24709, 37267,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 54613, 0, 3,
                                                                       54193, 36889, 54221,
                                                                       24709, 24754, 37330,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 54697, 0, 3,
                                                                       54221, 36910, 54249,
                                                                       24754, 24799, 37393,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 54781, 0, 3,
                                                                       54249, 36931, 54277,
                                                                       24799, 24844, 37456,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 54865, 0, 3,
                                                                       54277, 36952, 54305,
                                                                       24844, 24889, 37519,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 54949, 0, 3,
                                                                       54305, 36973, 54333,
                                                                       24889, 24934, 37582,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 55033, 0, 3,
                                                                       54361, 37141, 54445,
                                                                       25024, 25114, 37897,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 55201, 0, 3,
                                                                       54445, 37204, 54529,
                                                                       25114, 25204, 38023,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 55369, 0, 3,
                                                                       54529, 37267, 54613,
                                                                       25204, 25294, 38149,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 55537, 0, 3,
                                                                       54613, 37330, 54697,
                                                                       25294, 25384, 38275,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 55705, 0, 3,
                                                                       54697, 37393, 54781,
                                                                       25384, 25474, 38401,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 55873, 0, 3,
                                                                       54781, 37456, 54865,
                                                                       25474, 25564, 38527,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 56041, 0, 3,
                                                                       54865, 37519, 54949,
                                                                       25564, 25654, 38653,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 56209, 0, 3,
                                                                       55033, 37897, 55201,
                                                                       25834, 25984, 39199,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 56489, 0, 3,
                                                                       55201, 38023, 55369,
                                                                       25984, 26134, 39409,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 56769, 0, 3,
                                                                       55369, 38149, 55537,
                                                                       26134, 26284, 39619,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 57049, 0, 3,
                                                                       55537, 38275, 55705,
                                                                       26284, 26434, 39829,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 57329, 0, 3,
                                                                       55705, 38401, 55873,
                                                                       26434, 26584, 40039,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 57609, 0, 3,
                                                                       55873, 38527, 56041,
                                                                       26584, 26734, 40249,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 57889, 0, 3,
                                                                       56209, 39199, 56489,
                                                                       27034, 27259, 41089,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 58309, 0, 3,
                                                                       56489, 39409, 56769,
                                                                       27259, 27484, 41404,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 58729, 0, 3,
                                                                       56769, 39619, 57049,
                                                                       27484, 27709, 41719,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 59149, 0, 3,
                                                                       57049, 39829, 57329,
                                                                       27709, 27934, 42034,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 59569, 0, 3,
                                                                       57329, 40039, 57609,
                                                                       27934, 28159, 42349,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 59989, 0, 3,
                                                                       57889, 41089, 58309,
                                                                       28609, 28924, 43546,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 60577, 0, 3,
                                                                       58309, 41404, 58729,
                                                                       28924, 29239, 43987,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 61165, 0, 3,
                                                                       58729, 41719, 59149,
                                                                       29239, 29554, 44428,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 61753, 0, 3,
                                                                       59149, 42034, 59569,
                                                                       29554, 29869, 44869,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 62341, 0, 3,
                                                                       59989, 43546, 60577,
                                                                       30499, 30919, 46486,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 63125, 0, 3,
                                                                       60577, 43987, 61165,
                                                                       30919, 31339, 47074,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 63909, 0, 3,
                                                                       61165, 44428, 61753,
                                                                       31339, 31759, 47662,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 64693, 0, 3,
                                                                       62341, 46486, 63125,
                                                                       32599, 33139, 49762,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 65701, 0, 3,
                                                                       63125, 47074, 63909,
                                                                       33139, 33679, 50518,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 66709, 0, 3,
                                                                       64693, 49762, 65701,
                                                                       34759, 35434, 53164,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 67969, 3, 36784,
                                                                       36805, 54109, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 68005, 3, 36805,
                                                                       36826, 54137, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 68041, 3, 36826,
                                                                       36847, 54165, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 68077, 3, 36847,
                                                                       36868, 54193, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 68113, 3, 36868,
                                                                       36889, 54221, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 68149, 3, 36889,
                                                                       36910, 54249, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 68185, 3, 36910,
                                                                       36931, 54277, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 68221, 3, 36931,
                                                                       36952, 54305, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 68257, 3, 36952,
                                                                       36973, 54333, ncols,
                                                                       gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 68293, 0, 3,
                                                                       67969, 54109, 68005,
                                                                       37015, 37078, 54361,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 68401, 0, 3,
                                                                       68005, 54137, 68041,
                                                                       37078, 37141, 54445,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 68509, 0, 3,
                                                                       68041, 54165, 68077,
                                                                       37141, 37204, 54529,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 68617, 0, 3,
                                                                       68077, 54193, 68113,
                                                                       37204, 37267, 54613,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 68725, 0, 3,
                                                                       68113, 54221, 68149,
                                                                       37267, 37330, 54697,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 68833, 0, 3,
                                                                       68149, 54249, 68185,
                                                                       37330, 37393, 54781,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 68941, 0, 3,
                                                                       68185, 54277, 68221,
                                                                       37393, 37456, 54865,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 69049, 0, 3,
                                                                       68221, 54305, 68257,
                                                                       37456, 37519, 54949,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 69157, 0, 3,
                                                                       68293, 54361, 68401,
                                                                       37645, 37771, 55033,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 69373, 0, 3,
                                                                       68401, 54445, 68509,
                                                                       37771, 37897, 55201,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 69589, 0, 3,
                                                                       68509, 54529, 68617,
                                                                       37897, 38023, 55369,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 69805, 0, 3,
                                                                       68617, 54613, 68725,
                                                                       38023, 38149, 55537,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 70021, 0, 3,
                                                                       68725, 54697, 68833,
                                                                       38149, 38275, 55705,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 70237, 0, 3,
                                                                       68833, 54781, 68941,
                                                                       38275, 38401, 55873,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 70453, 0, 3,
                                                                       68941, 54865, 69049,
                                                                       38401, 38527, 56041,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 70669, 0, 3,
                                                                       69157, 55033, 69373,
                                                                       38779, 38989, 56209,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 71029, 0, 3,
                                                                       69373, 55201, 69589,
                                                                       38989, 39199, 56489,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 71389, 0, 3,
                                                                       69589, 55369, 69805,
                                                                       39199, 39409, 56769,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 71749, 0, 3,
                                                                       69805, 55537, 70021,
                                                                       39409, 39619, 57049,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 72109, 0, 3,
                                                                       70021, 55705, 70237,
                                                                       39619, 39829, 57329,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 72469, 0, 3,
                                                                       70237, 55873, 70453,
                                                                       39829, 40039, 57609,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 72829, 0, 3,
                                                                       70669, 56209, 71029,
                                                                       40459, 40774, 57889,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 73369, 0, 3,
                                                                       71029, 56489, 71389,
                                                                       40774, 41089, 58309,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 73909, 0, 3,
                                                                       71389, 56769, 71749,
                                                                       41089, 41404, 58729,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 74449, 0, 3,
                                                                       71749, 57049, 72109,
                                                                       41404, 41719, 59149,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 74989, 0, 3,
                                                                       72109, 57329, 72469,
                                                                       41719, 42034, 59569,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 75529, 0, 3,
                                                                       72829, 57889, 73369,
                                                                       42664, 43105, 59989,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 76285, 0, 3,
                                                                       73369, 58309, 73909,
                                                                       43105, 43546, 60577,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 77041, 0, 3,
                                                                       73909, 58729, 74449,
                                                                       43546, 43987, 61165,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 77797, 0, 3,
                                                                       74449, 59149, 74989,
                                                                       43987, 44428, 61753,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 78553, 0, 3,
                                                                       75529, 59989, 76285,
                                                                       45310, 45898, 62341,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 79561, 0, 3,
                                                                       76285, 60577, 77041,
                                                                       45898, 46486, 63125,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 80569, 0, 3,
                                                                       77041, 61165, 77797,
                                                                       46486, 47074, 63909,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 81577, 0, 3,
                                                                       78553, 62341, 79561,
                                                                       48250, 49006, 64693,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 82873, 0, 3,
                                                                       79561, 63125, 80569,
                                                                       49006, 49762, 65701,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 84169, 0, 3,
                                                                       81577, 64693, 82873,
                                                                       51274, 52219, 66709,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 85789, 75529, 756, ncols);

                    simdfunc::contract_primitives(buffer, 86860, 78553, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 88288, 81577, 1296, ncols);

                    simdfunc::contract_primitives(buffer, 90124, 84169, 1620, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 86545, 85789, 21, 1, nmax);

        simdtrf::transform_k_inner(buffer, 87868, 86860, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 89584, 88288, 36, 1, nmax);

        simdtrf::transform_k_inner(buffer, 91744, 90124, 45, 1, nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 92419, 86545, 87868, 15, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 93364, 87868, 89584, 15, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 94624, 89584, 91744, 15, nmax);

        simdtrf::compute_hrr_dh(buffer, coordinates, 96244, 92419, 93364, 15, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 98134, 93364, 94624, 15, nmax);

        simdtrf::compute_hrr_fh(buffer, coordinates, 100654, 96244, 98134, 15, nmax);

        simdtrf::transform_h_inner(buffer, 103804, 100654, 10, 15, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 103804, 165, nmax);
    }

    for (size_t m = 0; m < 1155; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
