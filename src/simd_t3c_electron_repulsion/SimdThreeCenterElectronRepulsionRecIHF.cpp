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


#include "SimdThreeCenterElectronRepulsionRecIHF.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

#include "SimdThreeCenterElectronRepulsionVrrRecDSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferID.hpp"
#include "SimdTransferIF.hpp"
#include "SimdTransferIG.hpp"
#include "SimdTransferIH.hpp"
#include "SimdTransferIP.hpp"
#include "SimdTransferKD.hpp"
#include "SimdTransferKF.hpp"
#include "SimdTransferKG.hpp"
#include "SimdTransferKP.hpp"
#include "SimdTransferLD.hpp"
#include "SimdTransferLF.hpp"
#include "SimdTransferLP.hpp"
#include "SimdTransferMD.hpp"
#include "SimdTransferMP.hpp"
#include "SimdTransferNP.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_ihf_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_ihf_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 67795, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1001 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 67795, 30219, 4690, dimensions);

        for (size_t i = 0; i < nprim_a; i++)
        {
            for (size_t j = 0; j < nprim_b; j++)
            {
                const auto p = a_exps[i] + b_exps[j];

                const auto mu = a_exps[i] * b_exps[j] / p;

                const auto fovl = a_norms[i] * b_norms[j];

                const auto fa = -b_exps[j] / p;

                const auto fc = b_exps[j] / p;

                simdfunc::compute_pa(buffer, coordinates, 0, nmax, fa);

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
                                                        5, 6, 7, 8, 9, 10, 11, 12, 13, 14},
                                                        ncols, fj, mu, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 21, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 24, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 27, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 30, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 33, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 36, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 39, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 42, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 45, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 48, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 51, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 54, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 57, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 60, 0, 3, 7, 8,
                                                                       21, 24, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 66, 0, 3, 8, 9,
                                                                       24, 27, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 72, 0, 3, 9, 10,
                                                                       27, 30, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 78, 0, 3, 10, 11,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 84, 0, 3, 11, 12,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 90, 0, 3, 12, 13,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 96, 0, 3, 13, 14,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 102, 0, 3, 14, 15,
                                                                       42, 45, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 108, 0, 3, 15, 16,
                                                                       45, 48, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 114, 0, 3, 16, 17,
                                                                       48, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 120, 0, 3, 17, 18,
                                                                       51, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 126, 0, 3, 18, 19,
                                                                       54, 57, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 132, 0, 3, 21, 24,
                                                                       60, 66, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 142, 0, 3, 24, 27,
                                                                       66, 72, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 152, 0, 3, 27, 30,
                                                                       72, 78, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 162, 0, 3, 30, 33,
                                                                       78, 84, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 172, 0, 3, 33, 36,
                                                                       84, 90, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 182, 0, 3, 36, 39,
                                                                       90, 96, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 192, 0, 3, 39, 42,
                                                                       96, 102, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 202, 0, 3, 42, 45,
                                                                       102, 108, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 212, 0, 3, 45, 48,
                                                                       108, 114, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 222, 0, 3, 48, 51,
                                                                       114, 120, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 232, 0, 3, 51, 54,
                                                                       120, 126, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 242, 0, 3, 60, 66,
                                                                       132, 142, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 257, 0, 3, 66, 72,
                                                                       142, 152, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 272, 0, 3, 72, 78,
                                                                       152, 162, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 287, 0, 3, 78, 84,
                                                                       162, 172, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 302, 0, 3, 84, 90,
                                                                       172, 182, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 317, 0, 3, 90, 96,
                                                                       182, 192, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 332, 0, 3, 96,
                                                                       102, 192, 202, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 347, 0, 3, 102,
                                                                       108, 202, 212, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 362, 0, 3, 108,
                                                                       114, 212, 222, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 377, 0, 3, 114,
                                                                       120, 222, 232, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 392, 0, 3, 132,
                                                                       142, 242, 257, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 413, 0, 3, 142,
                                                                       152, 257, 272, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 434, 0, 3, 152,
                                                                       162, 272, 287, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 455, 0, 3, 162,
                                                                       172, 287, 302, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 476, 0, 3, 172,
                                                                       182, 302, 317, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 497, 0, 3, 182,
                                                                       192, 317, 332, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 518, 0, 3, 192,
                                                                       202, 332, 347, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 539, 0, 3, 202,
                                                                       212, 347, 362, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 560, 0, 3, 212,
                                                                       222, 362, 377, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 581, 0, 3, 242,
                                                                       257, 392, 413, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 609, 0, 3, 257,
                                                                       272, 413, 434, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 637, 0, 3, 272,
                                                                       287, 434, 455, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 665, 0, 3, 287,
                                                                       302, 455, 476, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 693, 0, 3, 302,
                                                                       317, 476, 497, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 721, 0, 3, 317,
                                                                       332, 497, 518, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 749, 0, 3, 332,
                                                                       347, 518, 539, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 777, 0, 3, 347,
                                                                       362, 539, 560, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 805, 0, 3, 392,
                                                                       413, 581, 609, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 841, 0, 3, 413,
                                                                       434, 609, 637, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 877, 0, 3, 434,
                                                                       455, 637, 665, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 913, 0, 3, 455,
                                                                       476, 665, 693, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 949, 0, 3, 476,
                                                                       497, 693, 721, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 985, 0, 3, 497,
                                                                       518, 721, 749, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1021, 0, 3, 518,
                                                                       539, 749, 777, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1057, 0, 3, 581,
                                                                       609, 805, 841, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1102, 0, 3, 609,
                                                                       637, 841, 877, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1147, 0, 3, 637,
                                                                       665, 877, 913, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1192, 0, 3, 665,
                                                                       693, 913, 949, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1237, 0, 3, 693,
                                                                       721, 949, 985, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1282, 0, 3, 721,
                                                                       749, 985, 1021, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1327, 0, 3, 805,
                                                                       841, 1057, 1102, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1382, 0, 3, 841,
                                                                       877, 1102, 1147, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1437, 0, 3, 877,
                                                                       913, 1147, 1192, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1492, 0, 3, 913,
                                                                       949, 1192, 1237, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1547, 0, 3, 949,
                                                                       985, 1237, 1282, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 1602, 0, 3, 1057,
                                                                       1102, 1327, 1382, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 1668, 0, 3, 1102,
                                                                       1147, 1382, 1437, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 1734, 0, 3, 1147,
                                                                       1192, 1437, 1492, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 1800, 0, 3, 1192,
                                                                       1237, 1492, 1547, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 1866, 0, 3, 1327,
                                                                       1382, 1602, 1668, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 1944, 0, 3, 1382,
                                                                       1437, 1668, 1734, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 2022, 0, 3, 1437,
                                                                       1492, 1734, 1800, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2100, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2103, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2106, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2109, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2112, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2115, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2118, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2121, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2124, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2127, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2130, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2133, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2136, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2139, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2142, 3, 7, 21,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2151, 3, 8, 24,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2160, 3, 9, 27,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2169, 3, 10, 30,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2178, 3, 11, 33,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2187, 3, 12, 36,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2196, 3, 13, 39,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2205, 3, 14, 42,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2214, 3, 15, 45,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2223, 3, 16, 48,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2232, 3, 17, 51,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2241, 3, 18, 54,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2250, 3, 19, 57,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2259, 3, 21, 60,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2277, 3, 24, 66,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2295, 3, 27, 72,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2313, 3, 30, 78,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2331, 3, 33, 84,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2349, 3, 36, 90,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2367, 3, 39, 96,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2385, 3, 42, 102,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2403, 3, 45, 108,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2421, 3, 48, 114,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2439, 3, 51, 120,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2457, 3, 54, 126,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2475, 3, 60, 132,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2505, 3, 66, 142,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2535, 3, 72, 152,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2565, 3, 78, 162,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2595, 3, 84, 172,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2625, 3, 90, 182,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2655, 3, 96, 192,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2685, 3, 102, 202,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2715, 3, 108, 212,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2745, 3, 114, 222,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2775, 3, 120, 232,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2805, 3, 132, 242,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2850, 3, 142, 257,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2895, 3, 152, 272,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2940, 3, 162, 287,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2985, 3, 172, 302,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3030, 3, 182, 317,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3075, 3, 192, 332,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3120, 3, 202, 347,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3165, 3, 212, 362,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3210, 3, 222, 377,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3255, 3, 242, 392,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3318, 3, 257, 413,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3381, 3, 272, 434,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3444, 3, 287, 455,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3507, 3, 302, 476,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3570, 3, 317, 497,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3633, 3, 332, 518,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3696, 3, 347, 539,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3759, 3, 362, 560,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3822, 3, 392, 581,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3906, 3, 413, 609,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3990, 3, 434, 637,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4074, 3, 455, 665,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4158, 3, 476, 693,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4242, 3, 497, 721,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4326, 3, 518, 749,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4410, 3, 539, 777,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4494, 3, 581, 805,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4602, 3, 609, 841,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4710, 3, 637, 877,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4818, 3, 665, 913,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4926, 3, 693, 949,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5034, 3, 721, 985,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5142, 3, 749,
                                                                       1021, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5250, 3, 805,
                                                                       1057, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5385, 3, 841,
                                                                       1102, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5520, 3, 877,
                                                                       1147, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5655, 3, 913,
                                                                       1192, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5790, 3, 949,
                                                                       1237, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5925, 3, 985,
                                                                       1282, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 6060, 3, 1057,
                                                                       1327, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 6225, 3, 1102,
                                                                       1382, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 6390, 3, 1147,
                                                                       1437, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 6555, 3, 1192,
                                                                       1492, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 6720, 3, 1237,
                                                                       1547, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 6885, 3, 1327,
                                                                       1602, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 7083, 3, 1382,
                                                                       1668, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 7281, 3, 1437,
                                                                       1734, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 7479, 3, 1492,
                                                                       1800, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 7677, 3, 1602,
                                                                       1866, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 7911, 3, 1668,
                                                                       1944, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 8145, 3, 1734,
                                                                       2022, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8379, 3, 7, 8,
                                                                       2106, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8385, 3, 8, 9,
                                                                       2109, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8391, 3, 9, 10,
                                                                       2112, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8397, 3, 10, 11,
                                                                       2115, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8403, 3, 11, 12,
                                                                       2118, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8409, 3, 12, 13,
                                                                       2121, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8415, 3, 13, 14,
                                                                       2124, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8421, 3, 14, 15,
                                                                       2127, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8427, 3, 15, 16,
                                                                       2130, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8433, 3, 16, 17,
                                                                       2133, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8439, 3, 17, 18,
                                                                       2136, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8445, 3, 18, 19,
                                                                       2139, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8451, 0, 3, 8379,
                                                                       2106, 8385, 2160, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8469, 0, 3, 8385,
                                                                       2109, 8391, 2169, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8487, 0, 3, 8391,
                                                                       2112, 8397, 2178, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8505, 0, 3, 8397,
                                                                       2115, 8403, 2187, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8523, 0, 3, 8403,
                                                                       2118, 8409, 2196, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8541, 0, 3, 8409,
                                                                       2121, 8415, 2205, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8559, 0, 3, 8415,
                                                                       2124, 8421, 2214, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8577, 0, 3, 8421,
                                                                       2127, 8427, 2223, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8595, 0, 3, 8427,
                                                                       2130, 8433, 2232, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8613, 0, 3, 8433,
                                                                       2133, 8439, 2241, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8631, 0, 3, 8439,
                                                                       2136, 8445, 2250, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8649, 0, 3, 8451,
                                                                       2160, 8469, 60, 66, 2295,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8685, 0, 3, 8469,
                                                                       2169, 8487, 66, 72, 2313,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8721, 0, 3, 8487,
                                                                       2178, 8505, 72, 78, 2331,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8757, 0, 3, 8505,
                                                                       2187, 8523, 78, 84, 2349,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8793, 0, 3, 8523,
                                                                       2196, 8541, 84, 90, 2367,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8829, 0, 3, 8541,
                                                                       2205, 8559, 90, 96, 2385,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8865, 0, 3, 8559,
                                                                       2214, 8577, 96, 102, 2403,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8901, 0, 3, 8577,
                                                                       2223, 8595, 102, 108,
                                                                       2421, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8937, 0, 3, 8595,
                                                                       2232, 8613, 108, 114,
                                                                       2439, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8973, 0, 3, 8613,
                                                                       2241, 8631, 114, 120,
                                                                       2457, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9009, 0, 3, 8649,
                                                                       2295, 8685, 132, 142,
                                                                       2535, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9069, 0, 3, 8685,
                                                                       2313, 8721, 142, 152,
                                                                       2565, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9129, 0, 3, 8721,
                                                                       2331, 8757, 152, 162,
                                                                       2595, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9189, 0, 3, 8757,
                                                                       2349, 8793, 162, 172,
                                                                       2625, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9249, 0, 3, 8793,
                                                                       2367, 8829, 172, 182,
                                                                       2655, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9309, 0, 3, 8829,
                                                                       2385, 8865, 182, 192,
                                                                       2685, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9369, 0, 3, 8865,
                                                                       2403, 8901, 192, 202,
                                                                       2715, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9429, 0, 3, 8901,
                                                                       2421, 8937, 202, 212,
                                                                       2745, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9489, 0, 3, 8937,
                                                                       2439, 8973, 212, 222,
                                                                       2775, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 9549, 0, 3, 9009,
                                                                       2535, 9069, 242, 257,
                                                                       2895, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 9639, 0, 3, 9069,
                                                                       2565, 9129, 257, 272,
                                                                       2940, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 9729, 0, 3, 9129,
                                                                       2595, 9189, 272, 287,
                                                                       2985, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 9819, 0, 3, 9189,
                                                                       2625, 9249, 287, 302,
                                                                       3030, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 9909, 0, 3, 9249,
                                                                       2655, 9309, 302, 317,
                                                                       3075, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 9999, 0, 3, 9309,
                                                                       2685, 9369, 317, 332,
                                                                       3120, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10089, 0, 3, 9369,
                                                                       2715, 9429, 332, 347,
                                                                       3165, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10179, 0, 3, 9429,
                                                                       2745, 9489, 347, 362,
                                                                       3210, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 10269, 0, 3, 9549,
                                                                       2895, 9639, 392, 413,
                                                                       3381, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 10395, 0, 3, 9639,
                                                                       2940, 9729, 413, 434,
                                                                       3444, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 10521, 0, 3, 9729,
                                                                       2985, 9819, 434, 455,
                                                                       3507, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 10647, 0, 3, 9819,
                                                                       3030, 9909, 455, 476,
                                                                       3570, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 10773, 0, 3, 9909,
                                                                       3075, 9999, 476, 497,
                                                                       3633, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 10899, 0, 3, 9999,
                                                                       3120, 10089, 497, 518,
                                                                       3696, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 11025, 0, 3,
                                                                       10089, 3165, 10179, 518,
                                                                       539, 3759, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 11151, 0, 3,
                                                                       10269, 3381, 10395, 581,
                                                                       609, 3990, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 11319, 0, 3,
                                                                       10395, 3444, 10521, 609,
                                                                       637, 4074, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 11487, 0, 3,
                                                                       10521, 3507, 10647, 637,
                                                                       665, 4158, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 11655, 0, 3,
                                                                       10647, 3570, 10773, 665,
                                                                       693, 4242, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 11823, 0, 3,
                                                                       10773, 3633, 10899, 693,
                                                                       721, 4326, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 11991, 0, 3,
                                                                       10899, 3696, 11025, 721,
                                                                       749, 4410, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 12159, 0, 3,
                                                                       11151, 3990, 11319, 805,
                                                                       841, 4710, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 12375, 0, 3,
                                                                       11319, 4074, 11487, 841,
                                                                       877, 4818, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 12591, 0, 3,
                                                                       11487, 4158, 11655, 877,
                                                                       913, 4926, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 12807, 0, 3,
                                                                       11655, 4242, 11823, 913,
                                                                       949, 5034, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 13023, 0, 3,
                                                                       11823, 4326, 11991, 949,
                                                                       985, 5142, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 13239, 0, 3,
                                                                       12159, 4710, 12375, 1057,
                                                                       1102, 5520, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 13509, 0, 3,
                                                                       12375, 4818, 12591, 1102,
                                                                       1147, 5655, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 13779, 0, 3,
                                                                       12591, 4926, 12807, 1147,
                                                                       1192, 5790, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 14049, 0, 3,
                                                                       12807, 5034, 13023, 1192,
                                                                       1237, 5925, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 14319, 0, 3,
                                                                       13239, 5520, 13509, 1327,
                                                                       1382, 6390, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 14649, 0, 3,
                                                                       13509, 5655, 13779, 1382,
                                                                       1437, 6555, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 14979, 0, 3,
                                                                       13779, 5790, 14049, 1437,
                                                                       1492, 6720, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 15309, 0, 3,
                                                                       14319, 6390, 14649, 1602,
                                                                       1668, 7281, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 15705, 0, 3,
                                                                       14649, 6555, 14979, 1668,
                                                                       1734, 7479, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 16101, 0, 3,
                                                                       15309, 7281, 15705, 1866,
                                                                       1944, 8145, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16569, 3, 2100,
                                                                       2103, 8379, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16579, 3, 2103,
                                                                       2106, 8385, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16589, 3, 2106,
                                                                       2109, 8391, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16599, 3, 2109,
                                                                       2112, 8397, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16609, 3, 2112,
                                                                       2115, 8403, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16619, 3, 2115,
                                                                       2118, 8409, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16629, 3, 2118,
                                                                       2121, 8415, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16639, 3, 2121,
                                                                       2124, 8421, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16649, 3, 2124,
                                                                       2127, 8427, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16659, 3, 2127,
                                                                       2130, 8433, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16669, 3, 2130,
                                                                       2133, 8439, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16679, 3, 2133,
                                                                       2136, 8445, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 16689, 0, 3,
                                                                       16569, 8379, 16579, 2142,
                                                                       2151, 8451, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 16719, 0, 3,
                                                                       16579, 8385, 16589, 2151,
                                                                       2160, 8469, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 16749, 0, 3,
                                                                       16589, 8391, 16599, 2160,
                                                                       2169, 8487, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 16779, 0, 3,
                                                                       16599, 8397, 16609, 2169,
                                                                       2178, 8505, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 16809, 0, 3,
                                                                       16609, 8403, 16619, 2178,
                                                                       2187, 8523, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 16839, 0, 3,
                                                                       16619, 8409, 16629, 2187,
                                                                       2196, 8541, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 16869, 0, 3,
                                                                       16629, 8415, 16639, 2196,
                                                                       2205, 8559, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 16899, 0, 3,
                                                                       16639, 8421, 16649, 2205,
                                                                       2214, 8577, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 16929, 0, 3,
                                                                       16649, 8427, 16659, 2214,
                                                                       2223, 8595, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 16959, 0, 3,
                                                                       16659, 8433, 16669, 2223,
                                                                       2232, 8613, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 16989, 0, 3,
                                                                       16669, 8439, 16679, 2232,
                                                                       2241, 8631, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 17019, 0, 3,
                                                                       16689, 8451, 16719, 2259,
                                                                       2277, 8649, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 17079, 0, 3,
                                                                       16719, 8469, 16749, 2277,
                                                                       2295, 8685, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 17139, 0, 3,
                                                                       16749, 8487, 16779, 2295,
                                                                       2313, 8721, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 17199, 0, 3,
                                                                       16779, 8505, 16809, 2313,
                                                                       2331, 8757, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 17259, 0, 3,
                                                                       16809, 8523, 16839, 2331,
                                                                       2349, 8793, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 17319, 0, 3,
                                                                       16839, 8541, 16869, 2349,
                                                                       2367, 8829, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 17379, 0, 3,
                                                                       16869, 8559, 16899, 2367,
                                                                       2385, 8865, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 17439, 0, 3,
                                                                       16899, 8577, 16929, 2385,
                                                                       2403, 8901, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 17499, 0, 3,
                                                                       16929, 8595, 16959, 2403,
                                                                       2421, 8937, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 17559, 0, 3,
                                                                       16959, 8613, 16989, 2421,
                                                                       2439, 8973, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 17619, 0, 3,
                                                                       17019, 8649, 17079, 2475,
                                                                       2505, 9009, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 17719, 0, 3,
                                                                       17079, 8685, 17139, 2505,
                                                                       2535, 9069, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 17819, 0, 3,
                                                                       17139, 8721, 17199, 2535,
                                                                       2565, 9129, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 17919, 0, 3,
                                                                       17199, 8757, 17259, 2565,
                                                                       2595, 9189, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 18019, 0, 3,
                                                                       17259, 8793, 17319, 2595,
                                                                       2625, 9249, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 18119, 0, 3,
                                                                       17319, 8829, 17379, 2625,
                                                                       2655, 9309, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 18219, 0, 3,
                                                                       17379, 8865, 17439, 2655,
                                                                       2685, 9369, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 18319, 0, 3,
                                                                       17439, 8901, 17499, 2685,
                                                                       2715, 9429, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 18419, 0, 3,
                                                                       17499, 8937, 17559, 2715,
                                                                       2745, 9489, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 18519, 0, 3,
                                                                       17619, 9009, 17719, 2805,
                                                                       2850, 9549, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 18669, 0, 3,
                                                                       17719, 9069, 17819, 2850,
                                                                       2895, 9639, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 18819, 0, 3,
                                                                       17819, 9129, 17919, 2895,
                                                                       2940, 9729, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 18969, 0, 3,
                                                                       17919, 9189, 18019, 2940,
                                                                       2985, 9819, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 19119, 0, 3,
                                                                       18019, 9249, 18119, 2985,
                                                                       3030, 9909, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 19269, 0, 3,
                                                                       18119, 9309, 18219, 3030,
                                                                       3075, 9999, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 19419, 0, 3,
                                                                       18219, 9369, 18319, 3075,
                                                                       3120, 10089, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 19569, 0, 3,
                                                                       18319, 9429, 18419, 3120,
                                                                       3165, 10179, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 19719, 0, 3,
                                                                       18519, 9549, 18669, 3255,
                                                                       3318, 10269, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 19929, 0, 3,
                                                                       18669, 9639, 18819, 3318,
                                                                       3381, 10395, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 20139, 0, 3,
                                                                       18819, 9729, 18969, 3381,
                                                                       3444, 10521, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 20349, 0, 3,
                                                                       18969, 9819, 19119, 3444,
                                                                       3507, 10647, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 20559, 0, 3,
                                                                       19119, 9909, 19269, 3507,
                                                                       3570, 10773, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 20769, 0, 3,
                                                                       19269, 9999, 19419, 3570,
                                                                       3633, 10899, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 20979, 0, 3,
                                                                       19419, 10089, 19569, 3633,
                                                                       3696, 11025, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 21189, 0, 3,
                                                                       19719, 10269, 19929, 3822,
                                                                       3906, 11151, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 21469, 0, 3,
                                                                       19929, 10395, 20139, 3906,
                                                                       3990, 11319, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 21749, 0, 3,
                                                                       20139, 10521, 20349, 3990,
                                                                       4074, 11487, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 22029, 0, 3,
                                                                       20349, 10647, 20559, 4074,
                                                                       4158, 11655, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 22309, 0, 3,
                                                                       20559, 10773, 20769, 4158,
                                                                       4242, 11823, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 22589, 0, 3,
                                                                       20769, 10899, 20979, 4242,
                                                                       4326, 11991, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 22869, 0, 3,
                                                                       21189, 11151, 21469, 4494,
                                                                       4602, 12159, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 23229, 0, 3,
                                                                       21469, 11319, 21749, 4602,
                                                                       4710, 12375, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 23589, 0, 3,
                                                                       21749, 11487, 22029, 4710,
                                                                       4818, 12591, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 23949, 0, 3,
                                                                       22029, 11655, 22309, 4818,
                                                                       4926, 12807, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 24309, 0, 3,
                                                                       22309, 11823, 22589, 4926,
                                                                       5034, 13023, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 24669, 0, 3,
                                                                       22869, 12159, 23229, 5250,
                                                                       5385, 13239, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 25119, 0, 3,
                                                                       23229, 12375, 23589, 5385,
                                                                       5520, 13509, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 25569, 0, 3,
                                                                       23589, 12591, 23949, 5520,
                                                                       5655, 13779, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 26019, 0, 3,
                                                                       23949, 12807, 24309, 5655,
                                                                       5790, 14049, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 26469, 0, 3,
                                                                       24669, 13239, 25119, 6060,
                                                                       6225, 14319, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 27019, 0, 3,
                                                                       25119, 13509, 25569, 6225,
                                                                       6390, 14649, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 27569, 0, 3,
                                                                       25569, 13779, 26019, 6390,
                                                                       6555, 14979, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 28119, 0, 3,
                                                                       26469, 14319, 27019, 6885,
                                                                       7083, 15309, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 28779, 0, 3,
                                                                       27019, 14649, 27569, 7083,
                                                                       7281, 15705, ncols, gamma,
                                                                       p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 29439, 0, 3,
                                                                       28119, 15309, 28779, 7677,
                                                                       7911, 16101, ncols, gamma,
                                                                       p, q);

                    simdfunc::contract_primitives(buffer, 30219, 21189, 280, ncols);

                    simdfunc::contract_primitives(buffer, 30695, 22869, 360, ncols);

                    simdfunc::contract_primitives(buffer, 31307, 24669, 450, ncols);

                    simdfunc::contract_primitives(buffer, 32072, 26469, 550, ncols);

                    simdfunc::contract_primitives(buffer, 33007, 28119, 660, ncols);

                    simdfunc::contract_primitives(buffer, 34129, 29439, 780, ncols);
                }
            }
        }

        simdtrf::transform_f_inner(buffer, 30499, 30219, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 31055, 30695, 36, 1, nmax);

        simdtrf::transform_f_inner(buffer, 31757, 31307, 45, 1, nmax);

        simdtrf::transform_f_inner(buffer, 32622, 32072, 55, 1, nmax);

        simdtrf::transform_f_inner(buffer, 33667, 33007, 66, 1, nmax);

        simdtrf::transform_f_inner(buffer, 34909, 34129, 78, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 35455, 30499, 31055, 7, nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 36043, 31055, 31757, 7, nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 36799, 31757, 32622, 7, nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 37744, 32622, 33667, 7, nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 38899, 33667, 34909, 7, nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 40285, 35455, 36043, 7, nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 41461, 36043, 36799, 7, nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 42973, 36799, 37744, 7, nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 44863, 37744, 38899, 7, nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 47173, 40285, 41461, 7, nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 49133, 41461, 42973, 7, nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 51653, 42973, 44863, 7, nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 54803, 47173, 49133, 7, nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 57743, 49133, 51653, 7, nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 61523, 54803, 57743, 7, nmax);

        simdtrf::transform_h_inner(buffer, 65639, 61523, 28, 7, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 65639, 77, nmax);
    }

    for (size_t m = 0; m < 1001; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
