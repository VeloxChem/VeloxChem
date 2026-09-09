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


#include "SimdThreeCenterElectronRepulsionRecGIH.hpp"

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
compute_gih_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_gih_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 103544, 0, 0, dimensions);

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
        simdfunc::prepare_buffer(buffer, 103544, 69773, 6634, dimensions);

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

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1492, 0, 3, 889,
                                                                       925, 1177, 1222, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1547, 0, 3, 925,
                                                                       961, 1222, 1267, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1602, 0, 3, 961,
                                                                       997, 1267, 1312, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1657, 0, 3, 997,
                                                                       1033, 1312, 1357, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1712, 0, 3, 1033,
                                                                       1069, 1357, 1402, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1767, 0, 3, 1069,
                                                                       1105, 1402, 1447, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 1822, 0, 3, 1177,
                                                                       1222, 1492, 1547, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 1888, 0, 3, 1222,
                                                                       1267, 1547, 1602, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 1954, 0, 3, 1267,
                                                                       1312, 1602, 1657, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2020, 0, 3, 1312,
                                                                       1357, 1657, 1712, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2086, 0, 3, 1357,
                                                                       1402, 1712, 1767, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2152, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2155, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2158, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2161, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2164, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2167, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2170, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2173, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2176, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2179, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2182, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2185, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2188, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2191, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2194, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2197, 3, 9, 28,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2206, 3, 10, 31,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2215, 3, 11, 34,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2224, 3, 12, 37,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2233, 3, 13, 40,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2242, 3, 14, 43,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2251, 3, 15, 46,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2260, 3, 16, 49,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2269, 3, 17, 52,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2278, 3, 18, 55,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2287, 3, 19, 58,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2296, 3, 20, 61,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2305, 3, 22, 64,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2323, 3, 25, 70,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2341, 3, 28, 76,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2359, 3, 31, 82,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2377, 3, 34, 88,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2395, 3, 37, 94,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2413, 3, 40, 100,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2431, 3, 43, 106,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2449, 3, 46, 112,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2467, 3, 49, 118,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2485, 3, 52, 124,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2503, 3, 55, 130,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2521, 3, 58, 136,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2539, 3, 64, 142,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2569, 3, 70, 152,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2599, 3, 76, 162,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2629, 3, 82, 172,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2659, 3, 88, 182,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2689, 3, 94, 192,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2719, 3, 100, 202,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2749, 3, 106, 212,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2779, 3, 112, 222,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2809, 3, 118, 232,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2839, 3, 124, 242,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2869, 3, 130, 252,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2899, 3, 142, 262,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2944, 3, 152, 277,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2989, 3, 162, 292,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3034, 3, 172, 307,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3079, 3, 182, 322,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3124, 3, 192, 337,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3169, 3, 202, 352,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3214, 3, 212, 367,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3259, 3, 222, 382,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3304, 3, 232, 397,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3349, 3, 242, 412,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3394, 3, 262, 427,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3457, 3, 277, 448,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3520, 3, 292, 469,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3583, 3, 307, 490,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3646, 3, 322, 511,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3709, 3, 337, 532,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3772, 3, 352, 553,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3835, 3, 367, 574,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3898, 3, 382, 595,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3961, 3, 397, 616,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4024, 3, 427, 637,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4108, 3, 448, 665,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4192, 3, 469, 693,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4276, 3, 490, 721,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4360, 3, 511, 749,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4444, 3, 532, 777,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4528, 3, 553, 805,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4612, 3, 574, 833,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4696, 3, 595, 861,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4780, 3, 637, 889,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4888, 3, 665, 925,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4996, 3, 693, 961,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5104, 3, 721, 997,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5212, 3, 749,
                                                                       1033, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5320, 3, 777,
                                                                       1069, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5428, 3, 805,
                                                                       1105, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5536, 3, 833,
                                                                       1141, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5644, 3, 889,
                                                                       1177, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5779, 3, 925,
                                                                       1222, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5914, 3, 961,
                                                                       1267, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6049, 3, 997,
                                                                       1312, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6184, 3, 1033,
                                                                       1357, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6319, 3, 1069,
                                                                       1402, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6454, 3, 1105,
                                                                       1447, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 6589, 3, 1177,
                                                                       1492, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 6754, 3, 1222,
                                                                       1547, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 6919, 3, 1267,
                                                                       1602, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 7084, 3, 1312,
                                                                       1657, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 7249, 3, 1357,
                                                                       1712, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 7414, 3, 1402,
                                                                       1767, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 7579, 3, 1492,
                                                                       1822, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 7777, 3, 1547,
                                                                       1888, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 7975, 3, 1602,
                                                                       1954, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 8173, 3, 1657,
                                                                       2020, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 8371, 3, 1712,
                                                                       2086, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8569, 3, 7, 8,
                                                                       2158, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8575, 3, 8, 9,
                                                                       2161, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8581, 3, 9, 10,
                                                                       2164, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8587, 3, 10, 11,
                                                                       2167, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8593, 3, 11, 12,
                                                                       2170, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8599, 3, 12, 13,
                                                                       2173, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8605, 3, 13, 14,
                                                                       2176, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8611, 3, 14, 15,
                                                                       2179, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8617, 3, 15, 16,
                                                                       2182, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8623, 3, 16, 17,
                                                                       2185, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8629, 3, 17, 18,
                                                                       2188, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8635, 3, 18, 19,
                                                                       2191, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8641, 3, 19, 20,
                                                                       2194, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8647, 0, 3, 8569,
                                                                       2158, 8575, 2197, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8665, 0, 3, 8575,
                                                                       2161, 8581, 2206, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8683, 0, 3, 8581,
                                                                       2164, 8587, 2215, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8701, 0, 3, 8587,
                                                                       2167, 8593, 2224, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8719, 0, 3, 8593,
                                                                       2170, 8599, 2233, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8737, 0, 3, 8599,
                                                                       2173, 8605, 2242, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8755, 0, 3, 8605,
                                                                       2176, 8611, 2251, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8773, 0, 3, 8611,
                                                                       2179, 8617, 2260, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8791, 0, 3, 8617,
                                                                       2182, 8623, 2269, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8809, 0, 3, 8623,
                                                                       2185, 8629, 2278, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8827, 0, 3, 8629,
                                                                       2188, 8635, 2287, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8845, 0, 3, 8635,
                                                                       2191, 8641, 2296, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8863, 0, 3, 8647,
                                                                       2197, 8665, 64, 70, 2341,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8899, 0, 3, 8665,
                                                                       2206, 8683, 70, 76, 2359,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8935, 0, 3, 8683,
                                                                       2215, 8701, 76, 82, 2377,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8971, 0, 3, 8701,
                                                                       2224, 8719, 82, 88, 2395,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9007, 0, 3, 8719,
                                                                       2233, 8737, 88, 94, 2413,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9043, 0, 3, 8737,
                                                                       2242, 8755, 94, 100, 2431,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9079, 0, 3, 8755,
                                                                       2251, 8773, 100, 106,
                                                                       2449, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9115, 0, 3, 8773,
                                                                       2260, 8791, 106, 112,
                                                                       2467, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9151, 0, 3, 8791,
                                                                       2269, 8809, 112, 118,
                                                                       2485, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9187, 0, 3, 8809,
                                                                       2278, 8827, 118, 124,
                                                                       2503, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9223, 0, 3, 8827,
                                                                       2287, 8845, 124, 130,
                                                                       2521, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9259, 0, 3, 8863,
                                                                       2341, 8899, 142, 152,
                                                                       2599, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9319, 0, 3, 8899,
                                                                       2359, 8935, 152, 162,
                                                                       2629, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9379, 0, 3, 8935,
                                                                       2377, 8971, 162, 172,
                                                                       2659, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9439, 0, 3, 8971,
                                                                       2395, 9007, 172, 182,
                                                                       2689, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9499, 0, 3, 9007,
                                                                       2413, 9043, 182, 192,
                                                                       2719, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9559, 0, 3, 9043,
                                                                       2431, 9079, 192, 202,
                                                                       2749, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9619, 0, 3, 9079,
                                                                       2449, 9115, 202, 212,
                                                                       2779, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9679, 0, 3, 9115,
                                                                       2467, 9151, 212, 222,
                                                                       2809, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9739, 0, 3, 9151,
                                                                       2485, 9187, 222, 232,
                                                                       2839, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9799, 0, 3, 9187,
                                                                       2503, 9223, 232, 242,
                                                                       2869, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9859, 0, 3, 9259,
                                                                       2599, 9319, 262, 277,
                                                                       2989, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9949, 0, 3, 9319,
                                                                       2629, 9379, 277, 292,
                                                                       3034, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10039, 0, 3, 9379,
                                                                       2659, 9439, 292, 307,
                                                                       3079, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10129, 0, 3, 9439,
                                                                       2689, 9499, 307, 322,
                                                                       3124, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10219, 0, 3, 9499,
                                                                       2719, 9559, 322, 337,
                                                                       3169, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10309, 0, 3, 9559,
                                                                       2749, 9619, 337, 352,
                                                                       3214, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10399, 0, 3, 9619,
                                                                       2779, 9679, 352, 367,
                                                                       3259, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10489, 0, 3, 9679,
                                                                       2809, 9739, 367, 382,
                                                                       3304, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10579, 0, 3, 9739,
                                                                       2839, 9799, 382, 397,
                                                                       3349, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10669, 0, 3, 9859,
                                                                       2989, 9949, 427, 448,
                                                                       3520, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10795, 0, 3, 9949,
                                                                       3034, 10039, 448, 469,
                                                                       3583, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10921, 0, 3,
                                                                       10039, 3079, 10129, 469,
                                                                       490, 3646, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11047, 0, 3,
                                                                       10129, 3124, 10219, 490,
                                                                       511, 3709, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11173, 0, 3,
                                                                       10219, 3169, 10309, 511,
                                                                       532, 3772, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11299, 0, 3,
                                                                       10309, 3214, 10399, 532,
                                                                       553, 3835, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11425, 0, 3,
                                                                       10399, 3259, 10489, 553,
                                                                       574, 3898, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11551, 0, 3,
                                                                       10489, 3304, 10579, 574,
                                                                       595, 3961, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11677, 0, 3,
                                                                       10669, 3520, 10795, 637,
                                                                       665, 4192, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11845, 0, 3,
                                                                       10795, 3583, 10921, 665,
                                                                       693, 4276, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 12013, 0, 3,
                                                                       10921, 3646, 11047, 693,
                                                                       721, 4360, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 12181, 0, 3,
                                                                       11047, 3709, 11173, 721,
                                                                       749, 4444, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 12349, 0, 3,
                                                                       11173, 3772, 11299, 749,
                                                                       777, 4528, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 12517, 0, 3,
                                                                       11299, 3835, 11425, 777,
                                                                       805, 4612, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 12685, 0, 3,
                                                                       11425, 3898, 11551, 805,
                                                                       833, 4696, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 12853, 0, 3,
                                                                       11677, 4192, 11845, 889,
                                                                       925, 4996, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 13069, 0, 3,
                                                                       11845, 4276, 12013, 925,
                                                                       961, 5104, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 13285, 0, 3,
                                                                       12013, 4360, 12181, 961,
                                                                       997, 5212, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 13501, 0, 3,
                                                                       12181, 4444, 12349, 997,
                                                                       1033, 5320, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 13717, 0, 3,
                                                                       12349, 4528, 12517, 1033,
                                                                       1069, 5428, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 13933, 0, 3,
                                                                       12517, 4612, 12685, 1069,
                                                                       1105, 5536, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 14149, 0, 3,
                                                                       12853, 4996, 13069, 1177,
                                                                       1222, 5914, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 14419, 0, 3,
                                                                       13069, 5104, 13285, 1222,
                                                                       1267, 6049, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 14689, 0, 3,
                                                                       13285, 5212, 13501, 1267,
                                                                       1312, 6184, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 14959, 0, 3,
                                                                       13501, 5320, 13717, 1312,
                                                                       1357, 6319, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 15229, 0, 3,
                                                                       13717, 5428, 13933, 1357,
                                                                       1402, 6454, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 15499, 0, 3,
                                                                       14149, 5914, 14419, 1492,
                                                                       1547, 6919, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 15829, 0, 3,
                                                                       14419, 6049, 14689, 1547,
                                                                       1602, 7084, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 16159, 0, 3,
                                                                       14689, 6184, 14959, 1602,
                                                                       1657, 7249, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 16489, 0, 3,
                                                                       14959, 6319, 15229, 1657,
                                                                       1712, 7414, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 16819, 0, 3,
                                                                       15499, 6919, 15829, 1822,
                                                                       1888, 7975, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 17215, 0, 3,
                                                                       15829, 7084, 16159, 1888,
                                                                       1954, 8173, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 17611, 0, 3,
                                                                       16159, 7249, 16489, 1954,
                                                                       2020, 8371, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18007, 3, 2152,
                                                                       2155, 8569, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18017, 3, 2155,
                                                                       2158, 8575, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18027, 3, 2158,
                                                                       2161, 8581, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18037, 3, 2161,
                                                                       2164, 8587, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18047, 3, 2164,
                                                                       2167, 8593, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18057, 3, 2167,
                                                                       2170, 8599, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18067, 3, 2170,
                                                                       2173, 8605, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18077, 3, 2173,
                                                                       2176, 8611, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18087, 3, 2176,
                                                                       2179, 8617, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18097, 3, 2179,
                                                                       2182, 8623, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18107, 3, 2182,
                                                                       2185, 8629, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18117, 3, 2185,
                                                                       2188, 8635, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18127, 3, 2188,
                                                                       2191, 8641, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18137, 0, 3,
                                                                       18007, 8569, 18017, 8647,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18167, 0, 3,
                                                                       18017, 8575, 18027, 8665,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18197, 0, 3,
                                                                       18027, 8581, 18037, 8683,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18227, 0, 3,
                                                                       18037, 8587, 18047, 8701,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18257, 0, 3,
                                                                       18047, 8593, 18057, 8719,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18287, 0, 3,
                                                                       18057, 8599, 18067, 8737,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18317, 0, 3,
                                                                       18067, 8605, 18077, 8755,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18347, 0, 3,
                                                                       18077, 8611, 18087, 8773,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18377, 0, 3,
                                                                       18087, 8617, 18097, 8791,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18407, 0, 3,
                                                                       18097, 8623, 18107, 8809,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18437, 0, 3,
                                                                       18107, 8629, 18117, 8827,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18467, 0, 3,
                                                                       18117, 8635, 18127, 8845,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18497, 0, 3,
                                                                       18137, 8647, 18167, 2305,
                                                                       2323, 8863, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18557, 0, 3,
                                                                       18167, 8665, 18197, 2323,
                                                                       2341, 8899, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18617, 0, 3,
                                                                       18197, 8683, 18227, 2341,
                                                                       2359, 8935, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18677, 0, 3,
                                                                       18227, 8701, 18257, 2359,
                                                                       2377, 8971, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18737, 0, 3,
                                                                       18257, 8719, 18287, 2377,
                                                                       2395, 9007, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18797, 0, 3,
                                                                       18287, 8737, 18317, 2395,
                                                                       2413, 9043, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18857, 0, 3,
                                                                       18317, 8755, 18347, 2413,
                                                                       2431, 9079, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18917, 0, 3,
                                                                       18347, 8773, 18377, 2431,
                                                                       2449, 9115, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18977, 0, 3,
                                                                       18377, 8791, 18407, 2449,
                                                                       2467, 9151, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 19037, 0, 3,
                                                                       18407, 8809, 18437, 2467,
                                                                       2485, 9187, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 19097, 0, 3,
                                                                       18437, 8827, 18467, 2485,
                                                                       2503, 9223, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19157, 0, 3,
                                                                       18497, 8863, 18557, 2539,
                                                                       2569, 9259, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19257, 0, 3,
                                                                       18557, 8899, 18617, 2569,
                                                                       2599, 9319, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19357, 0, 3,
                                                                       18617, 8935, 18677, 2599,
                                                                       2629, 9379, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19457, 0, 3,
                                                                       18677, 8971, 18737, 2629,
                                                                       2659, 9439, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19557, 0, 3,
                                                                       18737, 9007, 18797, 2659,
                                                                       2689, 9499, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19657, 0, 3,
                                                                       18797, 9043, 18857, 2689,
                                                                       2719, 9559, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19757, 0, 3,
                                                                       18857, 9079, 18917, 2719,
                                                                       2749, 9619, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19857, 0, 3,
                                                                       18917, 9115, 18977, 2749,
                                                                       2779, 9679, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19957, 0, 3,
                                                                       18977, 9151, 19037, 2779,
                                                                       2809, 9739, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 20057, 0, 3,
                                                                       19037, 9187, 19097, 2809,
                                                                       2839, 9799, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 20157, 0, 3,
                                                                       19157, 9259, 19257, 2899,
                                                                       2944, 9859, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 20307, 0, 3,
                                                                       19257, 9319, 19357, 2944,
                                                                       2989, 9949, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 20457, 0, 3,
                                                                       19357, 9379, 19457, 2989,
                                                                       3034, 10039, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 20607, 0, 3,
                                                                       19457, 9439, 19557, 3034,
                                                                       3079, 10129, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 20757, 0, 3,
                                                                       19557, 9499, 19657, 3079,
                                                                       3124, 10219, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 20907, 0, 3,
                                                                       19657, 9559, 19757, 3124,
                                                                       3169, 10309, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 21057, 0, 3,
                                                                       19757, 9619, 19857, 3169,
                                                                       3214, 10399, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 21207, 0, 3,
                                                                       19857, 9679, 19957, 3214,
                                                                       3259, 10489, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 21357, 0, 3,
                                                                       19957, 9739, 20057, 3259,
                                                                       3304, 10579, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 21507, 0, 3,
                                                                       20157, 9859, 20307, 3394,
                                                                       3457, 10669, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 21717, 0, 3,
                                                                       20307, 9949, 20457, 3457,
                                                                       3520, 10795, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 21927, 0, 3,
                                                                       20457, 10039, 20607, 3520,
                                                                       3583, 10921, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 22137, 0, 3,
                                                                       20607, 10129, 20757, 3583,
                                                                       3646, 11047, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 22347, 0, 3,
                                                                       20757, 10219, 20907, 3646,
                                                                       3709, 11173, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 22557, 0, 3,
                                                                       20907, 10309, 21057, 3709,
                                                                       3772, 11299, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 22767, 0, 3,
                                                                       21057, 10399, 21207, 3772,
                                                                       3835, 11425, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 22977, 0, 3,
                                                                       21207, 10489, 21357, 3835,
                                                                       3898, 11551, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 23187, 0, 3,
                                                                       21507, 10669, 21717, 4024,
                                                                       4108, 11677, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 23467, 0, 3,
                                                                       21717, 10795, 21927, 4108,
                                                                       4192, 11845, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 23747, 0, 3,
                                                                       21927, 10921, 22137, 4192,
                                                                       4276, 12013, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 24027, 0, 3,
                                                                       22137, 11047, 22347, 4276,
                                                                       4360, 12181, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 24307, 0, 3,
                                                                       22347, 11173, 22557, 4360,
                                                                       4444, 12349, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 24587, 0, 3,
                                                                       22557, 11299, 22767, 4444,
                                                                       4528, 12517, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 24867, 0, 3,
                                                                       22767, 11425, 22977, 4528,
                                                                       4612, 12685, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 25147, 0, 3,
                                                                       23187, 11677, 23467, 4780,
                                                                       4888, 12853, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 25507, 0, 3,
                                                                       23467, 11845, 23747, 4888,
                                                                       4996, 13069, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 25867, 0, 3,
                                                                       23747, 12013, 24027, 4996,
                                                                       5104, 13285, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 26227, 0, 3,
                                                                       24027, 12181, 24307, 5104,
                                                                       5212, 13501, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 26587, 0, 3,
                                                                       24307, 12349, 24587, 5212,
                                                                       5320, 13717, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 26947, 0, 3,
                                                                       24587, 12517, 24867, 5320,
                                                                       5428, 13933, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 27307, 0, 3,
                                                                       25147, 12853, 25507, 5644,
                                                                       5779, 14149, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 27757, 0, 3,
                                                                       25507, 13069, 25867, 5779,
                                                                       5914, 14419, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 28207, 0, 3,
                                                                       25867, 13285, 26227, 5914,
                                                                       6049, 14689, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 28657, 0, 3,
                                                                       26227, 13501, 26587, 6049,
                                                                       6184, 14959, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 29107, 0, 3,
                                                                       26587, 13717, 26947, 6184,
                                                                       6319, 15229, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 29557, 0, 3,
                                                                       27307, 14149, 27757, 6589,
                                                                       6754, 15499, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 30107, 0, 3,
                                                                       27757, 14419, 28207, 6754,
                                                                       6919, 15829, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 30657, 0, 3,
                                                                       28207, 14689, 28657, 6919,
                                                                       7084, 16159, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 31207, 0, 3,
                                                                       28657, 14959, 29107, 7084,
                                                                       7249, 16489, ncols, gamma,
                                                                       p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 31757, 0, 3,
                                                                       29557, 15499, 30107, 7579,
                                                                       7777, 16819, ncols, gamma,
                                                                       p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 32417, 0, 3,
                                                                       30107, 15829, 30657, 7777,
                                                                       7975, 17215, ncols, gamma,
                                                                       p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 33077, 0, 3,
                                                                       30657, 16159, 31207, 7975,
                                                                       8173, 17611, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33737, 3, 8569,
                                                                       8575, 18027, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33752, 3, 8575,
                                                                       8581, 18037, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33767, 3, 8581,
                                                                       8587, 18047, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33782, 3, 8587,
                                                                       8593, 18057, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33797, 3, 8593,
                                                                       8599, 18067, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33812, 3, 8599,
                                                                       8605, 18077, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33827, 3, 8605,
                                                                       8611, 18087, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33842, 3, 8611,
                                                                       8617, 18097, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33857, 3, 8617,
                                                                       8623, 18107, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33872, 3, 8623,
                                                                       8629, 18117, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33887, 3, 8629,
                                                                       8635, 18127, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 33902, 0, 3,
                                                                       33737, 18027, 33752, 8647,
                                                                       8665, 18197, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 33947, 0, 3,
                                                                       33752, 18037, 33767, 8665,
                                                                       8683, 18227, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 33992, 0, 3,
                                                                       33767, 18047, 33782, 8683,
                                                                       8701, 18257, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 34037, 0, 3,
                                                                       33782, 18057, 33797, 8701,
                                                                       8719, 18287, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 34082, 0, 3,
                                                                       33797, 18067, 33812, 8719,
                                                                       8737, 18317, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 34127, 0, 3,
                                                                       33812, 18077, 33827, 8737,
                                                                       8755, 18347, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 34172, 0, 3,
                                                                       33827, 18087, 33842, 8755,
                                                                       8773, 18377, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 34217, 0, 3,
                                                                       33842, 18097, 33857, 8773,
                                                                       8791, 18407, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 34262, 0, 3,
                                                                       33857, 18107, 33872, 8791,
                                                                       8809, 18437, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 34307, 0, 3,
                                                                       33872, 18117, 33887, 8809,
                                                                       8827, 18467, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 34352, 0, 3,
                                                                       33902, 18197, 33947, 8863,
                                                                       8899, 18617, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 34442, 0, 3,
                                                                       33947, 18227, 33992, 8899,
                                                                       8935, 18677, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 34532, 0, 3,
                                                                       33992, 18257, 34037, 8935,
                                                                       8971, 18737, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 34622, 0, 3,
                                                                       34037, 18287, 34082, 8971,
                                                                       9007, 18797, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 34712, 0, 3,
                                                                       34082, 18317, 34127, 9007,
                                                                       9043, 18857, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 34802, 0, 3,
                                                                       34127, 18347, 34172, 9043,
                                                                       9079, 18917, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 34892, 0, 3,
                                                                       34172, 18377, 34217, 9079,
                                                                       9115, 18977, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 34982, 0, 3,
                                                                       34217, 18407, 34262, 9115,
                                                                       9151, 19037, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 35072, 0, 3,
                                                                       34262, 18437, 34307, 9151,
                                                                       9187, 19097, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 35162, 0, 3,
                                                                       34352, 18617, 34442, 9259,
                                                                       9319, 19357, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 35312, 0, 3,
                                                                       34442, 18677, 34532, 9319,
                                                                       9379, 19457, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 35462, 0, 3,
                                                                       34532, 18737, 34622, 9379,
                                                                       9439, 19557, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 35612, 0, 3,
                                                                       34622, 18797, 34712, 9439,
                                                                       9499, 19657, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 35762, 0, 3,
                                                                       34712, 18857, 34802, 9499,
                                                                       9559, 19757, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 35912, 0, 3,
                                                                       34802, 18917, 34892, 9559,
                                                                       9619, 19857, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 36062, 0, 3,
                                                                       34892, 18977, 34982, 9619,
                                                                       9679, 19957, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 36212, 0, 3,
                                                                       34982, 19037, 35072, 9679,
                                                                       9739, 20057, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 36362, 0, 3,
                                                                       35162, 19357, 35312, 9859,
                                                                       9949, 20457, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 36587, 0, 3,
                                                                       35312, 19457, 35462, 9949,
                                                                       10039, 20607, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 36812, 0, 3,
                                                                       35462, 19557, 35612,
                                                                       10039, 10129, 20757,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 37037, 0, 3,
                                                                       35612, 19657, 35762,
                                                                       10129, 10219, 20907,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 37262, 0, 3,
                                                                       35762, 19757, 35912,
                                                                       10219, 10309, 21057,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 37487, 0, 3,
                                                                       35912, 19857, 36062,
                                                                       10309, 10399, 21207,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 37712, 0, 3,
                                                                       36062, 19957, 36212,
                                                                       10399, 10489, 21357,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 37937, 0, 3,
                                                                       36362, 20457, 36587,
                                                                       10669, 10795, 21927,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 38252, 0, 3,
                                                                       36587, 20607, 36812,
                                                                       10795, 10921, 22137,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 38567, 0, 3,
                                                                       36812, 20757, 37037,
                                                                       10921, 11047, 22347,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 38882, 0, 3,
                                                                       37037, 20907, 37262,
                                                                       11047, 11173, 22557,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 39197, 0, 3,
                                                                       37262, 21057, 37487,
                                                                       11173, 11299, 22767,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 39512, 0, 3,
                                                                       37487, 21207, 37712,
                                                                       11299, 11425, 22977,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 39827, 0, 3,
                                                                       37937, 21927, 38252,
                                                                       11677, 11845, 23747,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 40247, 0, 3,
                                                                       38252, 22137, 38567,
                                                                       11845, 12013, 24027,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 40667, 0, 3,
                                                                       38567, 22347, 38882,
                                                                       12013, 12181, 24307,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 41087, 0, 3,
                                                                       38882, 22557, 39197,
                                                                       12181, 12349, 24587,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 41507, 0, 3,
                                                                       39197, 22767, 39512,
                                                                       12349, 12517, 24867,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 41927, 0, 3,
                                                                       39827, 23747, 40247,
                                                                       12853, 13069, 25867,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 42467, 0, 3,
                                                                       40247, 24027, 40667,
                                                                       13069, 13285, 26227,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 43007, 0, 3,
                                                                       40667, 24307, 41087,
                                                                       13285, 13501, 26587,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 43547, 0, 3,
                                                                       41087, 24587, 41507,
                                                                       13501, 13717, 26947,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 44087, 0, 3,
                                                                       41927, 25867, 42467,
                                                                       14149, 14419, 28207,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 44762, 0, 3,
                                                                       42467, 26227, 43007,
                                                                       14419, 14689, 28657,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 45437, 0, 3,
                                                                       43007, 26587, 43547,
                                                                       14689, 14959, 29107,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 46112, 0, 3,
                                                                       44087, 28207, 44762,
                                                                       15499, 15829, 30657,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 46937, 0, 3,
                                                                       44762, 28657, 45437,
                                                                       15829, 16159, 31207,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 47762, 0, 3,
                                                                       46112, 30657, 46937,
                                                                       16819, 17215, 33077,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 48752, 3, 18007,
                                                                       18017, 33737, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 48773, 3, 18017,
                                                                       18027, 33752, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 48794, 3, 18027,
                                                                       18037, 33767, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 48815, 3, 18037,
                                                                       18047, 33782, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 48836, 3, 18047,
                                                                       18057, 33797, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 48857, 3, 18057,
                                                                       18067, 33812, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 48878, 3, 18067,
                                                                       18077, 33827, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 48899, 3, 18077,
                                                                       18087, 33842, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 48920, 3, 18087,
                                                                       18097, 33857, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 48941, 3, 18097,
                                                                       18107, 33872, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 48962, 3, 18107,
                                                                       18117, 33887, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 48983, 0, 3,
                                                                       48752, 33737, 48773,
                                                                       18137, 18167, 33902,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 49046, 0, 3,
                                                                       48773, 33752, 48794,
                                                                       18167, 18197, 33947,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 49109, 0, 3,
                                                                       48794, 33767, 48815,
                                                                       18197, 18227, 33992,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 49172, 0, 3,
                                                                       48815, 33782, 48836,
                                                                       18227, 18257, 34037,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 49235, 0, 3,
                                                                       48836, 33797, 48857,
                                                                       18257, 18287, 34082,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 49298, 0, 3,
                                                                       48857, 33812, 48878,
                                                                       18287, 18317, 34127,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 49361, 0, 3,
                                                                       48878, 33827, 48899,
                                                                       18317, 18347, 34172,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 49424, 0, 3,
                                                                       48899, 33842, 48920,
                                                                       18347, 18377, 34217,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 49487, 0, 3,
                                                                       48920, 33857, 48941,
                                                                       18377, 18407, 34262,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 49550, 0, 3,
                                                                       48941, 33872, 48962,
                                                                       18407, 18437, 34307,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 49613, 0, 3,
                                                                       48983, 33902, 49046,
                                                                       18497, 18557, 34352,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 49739, 0, 3,
                                                                       49046, 33947, 49109,
                                                                       18557, 18617, 34442,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 49865, 0, 3,
                                                                       49109, 33992, 49172,
                                                                       18617, 18677, 34532,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 49991, 0, 3,
                                                                       49172, 34037, 49235,
                                                                       18677, 18737, 34622,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 50117, 0, 3,
                                                                       49235, 34082, 49298,
                                                                       18737, 18797, 34712,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 50243, 0, 3,
                                                                       49298, 34127, 49361,
                                                                       18797, 18857, 34802,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 50369, 0, 3,
                                                                       49361, 34172, 49424,
                                                                       18857, 18917, 34892,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 50495, 0, 3,
                                                                       49424, 34217, 49487,
                                                                       18917, 18977, 34982,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 50621, 0, 3,
                                                                       49487, 34262, 49550,
                                                                       18977, 19037, 35072,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 50747, 0, 3,
                                                                       49613, 34352, 49739,
                                                                       19157, 19257, 35162,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 50957, 0, 3,
                                                                       49739, 34442, 49865,
                                                                       19257, 19357, 35312,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 51167, 0, 3,
                                                                       49865, 34532, 49991,
                                                                       19357, 19457, 35462,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 51377, 0, 3,
                                                                       49991, 34622, 50117,
                                                                       19457, 19557, 35612,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 51587, 0, 3,
                                                                       50117, 34712, 50243,
                                                                       19557, 19657, 35762,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 51797, 0, 3,
                                                                       50243, 34802, 50369,
                                                                       19657, 19757, 35912,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 52007, 0, 3,
                                                                       50369, 34892, 50495,
                                                                       19757, 19857, 36062,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 52217, 0, 3,
                                                                       50495, 34982, 50621,
                                                                       19857, 19957, 36212,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 52427, 0, 3,
                                                                       50747, 35162, 50957,
                                                                       20157, 20307, 36362,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 52742, 0, 3,
                                                                       50957, 35312, 51167,
                                                                       20307, 20457, 36587,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 53057, 0, 3,
                                                                       51167, 35462, 51377,
                                                                       20457, 20607, 36812,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 53372, 0, 3,
                                                                       51377, 35612, 51587,
                                                                       20607, 20757, 37037,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 53687, 0, 3,
                                                                       51587, 35762, 51797,
                                                                       20757, 20907, 37262,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 54002, 0, 3,
                                                                       51797, 35912, 52007,
                                                                       20907, 21057, 37487,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 54317, 0, 3,
                                                                       52007, 36062, 52217,
                                                                       21057, 21207, 37712,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 54632, 0, 3,
                                                                       52427, 36362, 52742,
                                                                       21507, 21717, 37937,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 55073, 0, 3,
                                                                       52742, 36587, 53057,
                                                                       21717, 21927, 38252,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 55514, 0, 3,
                                                                       53057, 36812, 53372,
                                                                       21927, 22137, 38567,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 55955, 0, 3,
                                                                       53372, 37037, 53687,
                                                                       22137, 22347, 38882,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 56396, 0, 3,
                                                                       53687, 37262, 54002,
                                                                       22347, 22557, 39197,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 56837, 0, 3,
                                                                       54002, 37487, 54317,
                                                                       22557, 22767, 39512,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 57278, 0, 3,
                                                                       54632, 37937, 55073,
                                                                       23187, 23467, 39827,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 57866, 0, 3,
                                                                       55073, 38252, 55514,
                                                                       23467, 23747, 40247,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 58454, 0, 3,
                                                                       55514, 38567, 55955,
                                                                       23747, 24027, 40667,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 59042, 0, 3,
                                                                       55955, 38882, 56396,
                                                                       24027, 24307, 41087,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 59630, 0, 3,
                                                                       56396, 39197, 56837,
                                                                       24307, 24587, 41507,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 60218, 0, 3,
                                                                       57278, 39827, 57866,
                                                                       25147, 25507, 41927,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 60974, 0, 3,
                                                                       57866, 40247, 58454,
                                                                       25507, 25867, 42467,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 61730, 0, 3,
                                                                       58454, 40667, 59042,
                                                                       25867, 26227, 43007,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 62486, 0, 3,
                                                                       59042, 41087, 59630,
                                                                       26227, 26587, 43547,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 63242, 0, 3,
                                                                       60218, 41927, 60974,
                                                                       27307, 27757, 44087,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 64187, 0, 3,
                                                                       60974, 42467, 61730,
                                                                       27757, 28207, 44762,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 65132, 0, 3,
                                                                       61730, 43007, 62486,
                                                                       28207, 28657, 45437,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 66077, 0, 3,
                                                                       63242, 44087, 64187,
                                                                       29557, 30107, 46112,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 67232, 0, 3,
                                                                       64187, 44762, 65132,
                                                                       30107, 30657, 46937,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 68387, 0, 3,
                                                                       66077, 46112, 67232,
                                                                       31757, 32417, 47762,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 69773, 57278, 588, ncols);

                    simdfunc::contract_primitives(buffer, 70669, 60218, 756, ncols);

                    simdfunc::contract_primitives(buffer, 71821, 63242, 945, ncols);

                    simdfunc::contract_primitives(buffer, 73261, 66077, 1155, ncols);

                    simdfunc::contract_primitives(buffer, 75021, 68387, 1386, ncols);
                }
            }
        }

        simdtrf::transform_h_inner(buffer, 70361, 69773, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 71425, 70669, 36, 1, nmax);

        simdtrf::transform_h_inner(buffer, 72766, 71821, 45, 1, nmax);

        simdtrf::transform_h_inner(buffer, 74416, 73261, 55, 1, nmax);

        simdtrf::transform_h_inner(buffer, 76407, 75021, 66, 1, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 77133, 70361, 71425, 11, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 78057, 71425, 72766, 11, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 79245, 72766, 74416, 11, nmax);

        simdtrf::compute_hrr_pm(buffer, coordinates, 80730, 74416, 76407, 11, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 82545, 77133, 78057, 11, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 84393, 78057, 79245, 11, nmax);

        simdtrf::compute_hrr_dl(buffer, coordinates, 86769, 79245, 80730, 11, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 89739, 82545, 84393, 11, nmax);

        simdtrf::compute_hrr_fk(buffer, coordinates, 92819, 84393, 86769, 11, nmax);

        simdtrf::compute_hrr_gi(buffer, coordinates, 96779, 89739, 92819, 11, nmax);

        simdtrf::transform_i_inner(buffer, 101399, 96779, 15, 11, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 101399, 143, nmax);
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
