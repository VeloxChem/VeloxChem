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


#include "SimdThreeCenterElectronRepulsionRecHHD.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferHD.hpp"
#include "SimdTransferHF.hpp"
#include "SimdTransferHG.hpp"
#include "SimdTransferHH.hpp"
#include "SimdTransferHP.hpp"
#include "SimdTransferID.hpp"
#include "SimdTransferIF.hpp"
#include "SimdTransferIG.hpp"
#include "SimdTransferIP.hpp"
#include "SimdTransferKD.hpp"
#include "SimdTransferKF.hpp"
#include "SimdTransferKP.hpp"
#include "SimdTransferLD.hpp"
#include "SimdTransferLP.hpp"
#include "SimdTransferMP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformH.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_hhd_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_hhd_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 31311, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 605 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 31311, 10590, 2431, dimensions);

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

                simdfunc::compute_pair_exponent(buffer, coordinates, 6, nmax, mu);

                for (size_t k = 0; k < nprim_c; k++)
                {
                    const auto ncols = dimensions[(i * nprim_b + j) * nprim_c + k];

                    if (ncols == 0) continue;

                    const auto gamma = c_exps[k];

                    const auto q = p + gamma;

                    const auto fq = p * gamma / q;

                    const auto fj = 2.0 * fovl * c_norms[k] * pi * pi * std::sqrt(pi)
                                    / (p * gamma * std::sqrt(q));

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 7, 3, 12,
                                                             ncols, fj, 6, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 21, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 24, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 27, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 30, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 33, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 36, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 39, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 42, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 45, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 48, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 51, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 54, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 57, 0, 3, 8, 9,
                                                                       21, 24, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 63, 0, 3, 9, 10,
                                                                       24, 27, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 69, 0, 3, 10, 11,
                                                                       27, 30, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 75, 0, 3, 11, 12,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 81, 0, 3, 12, 13,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 87, 0, 3, 13, 14,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 93, 0, 3, 14, 15,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 99, 0, 3, 15, 16,
                                                                       42, 45, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 105, 0, 3, 16, 17,
                                                                       45, 48, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 111, 0, 3, 17, 18,
                                                                       48, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 117, 0, 3, 18, 19,
                                                                       51, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 123, 0, 3, 21, 24,
                                                                       57, 63, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 133, 0, 3, 24, 27,
                                                                       63, 69, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 143, 0, 3, 27, 30,
                                                                       69, 75, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 153, 0, 3, 30, 33,
                                                                       75, 81, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 163, 0, 3, 33, 36,
                                                                       81, 87, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 173, 0, 3, 36, 39,
                                                                       87, 93, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 183, 0, 3, 39, 42,
                                                                       93, 99, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 193, 0, 3, 42, 45,
                                                                       99, 105, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 203, 0, 3, 45, 48,
                                                                       105, 111, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 213, 0, 3, 48, 51,
                                                                       111, 117, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 223, 0, 3, 57, 63,
                                                                       123, 133, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 238, 0, 3, 63, 69,
                                                                       133, 143, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 253, 0, 3, 69, 75,
                                                                       143, 153, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 268, 0, 3, 75, 81,
                                                                       153, 163, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 283, 0, 3, 81, 87,
                                                                       163, 173, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 298, 0, 3, 87, 93,
                                                                       173, 183, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 313, 0, 3, 93, 99,
                                                                       183, 193, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 328, 0, 3, 99,
                                                                       105, 193, 203, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 343, 0, 3, 105,
                                                                       111, 203, 213, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 358, 0, 3, 123,
                                                                       133, 223, 238, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 379, 0, 3, 133,
                                                                       143, 238, 253, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 400, 0, 3, 143,
                                                                       153, 253, 268, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 421, 0, 3, 153,
                                                                       163, 268, 283, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 442, 0, 3, 163,
                                                                       173, 283, 298, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 463, 0, 3, 173,
                                                                       183, 298, 313, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 484, 0, 3, 183,
                                                                       193, 313, 328, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 505, 0, 3, 193,
                                                                       203, 328, 343, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 526, 0, 3, 223,
                                                                       238, 358, 379, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 554, 0, 3, 238,
                                                                       253, 379, 400, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 582, 0, 3, 253,
                                                                       268, 400, 421, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 610, 0, 3, 268,
                                                                       283, 421, 442, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 638, 0, 3, 283,
                                                                       298, 442, 463, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 666, 0, 3, 298,
                                                                       313, 463, 484, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 694, 0, 3, 313,
                                                                       328, 484, 505, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 722, 0, 3, 358,
                                                                       379, 526, 554, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 758, 0, 3, 379,
                                                                       400, 554, 582, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 794, 0, 3, 400,
                                                                       421, 582, 610, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 830, 0, 3, 421,
                                                                       442, 610, 638, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 866, 0, 3, 442,
                                                                       463, 638, 666, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 902, 0, 3, 463,
                                                                       484, 666, 694, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 938, 0, 3, 526,
                                                                       554, 722, 758, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 983, 0, 3, 554,
                                                                       582, 758, 794, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1028, 0, 3, 582,
                                                                       610, 794, 830, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1073, 0, 3, 610,
                                                                       638, 830, 866, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1118, 0, 3, 638,
                                                                       666, 866, 902, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1163, 0, 3, 722,
                                                                       758, 938, 983, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1218, 0, 3, 758,
                                                                       794, 983, 1028, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1273, 0, 3, 794,
                                                                       830, 1028, 1073, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1328, 0, 3, 830,
                                                                       866, 1073, 1118, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 1383, 0, 3, 938,
                                                                       983, 1163, 1218, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 1449, 0, 3, 983,
                                                                       1028, 1218, 1273, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 1515, 0, 3, 1028,
                                                                       1073, 1273, 1328, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1581, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1584, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1587, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1590, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1593, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1596, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1599, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1602, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1605, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1608, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1611, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1614, 3, 10, 27,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1623, 3, 11, 30,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1632, 3, 12, 33,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1641, 3, 13, 36,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1650, 3, 14, 39,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1659, 3, 15, 42,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1668, 3, 16, 45,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1677, 3, 17, 48,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1686, 3, 18, 51,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1695, 3, 19, 54,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1704, 3, 27, 69,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1722, 3, 30, 75,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1740, 3, 33, 81,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1758, 3, 36, 87,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1776, 3, 39, 93,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1794, 3, 42, 99,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1812, 3, 45, 105,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1830, 3, 48, 111,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1848, 3, 51, 117,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1866, 3, 69, 143,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1896, 3, 75, 153,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1926, 3, 81, 163,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1956, 3, 87, 173,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1986, 3, 93, 183,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2016, 3, 99, 193,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2046, 3, 105, 203,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2076, 3, 111, 213,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2106, 3, 143, 253,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2151, 3, 153, 268,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2196, 3, 163, 283,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2241, 3, 173, 298,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2286, 3, 183, 313,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2331, 3, 193, 328,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2376, 3, 203, 343,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2421, 3, 253, 400,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2484, 3, 268, 421,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2547, 3, 283, 442,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2610, 3, 298, 463,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2673, 3, 313, 484,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2736, 3, 328, 505,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2799, 3, 400, 582,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2883, 3, 421, 610,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2967, 3, 442, 638,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3051, 3, 463, 666,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3135, 3, 484, 694,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 3219, 3, 582, 794,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 3327, 3, 610, 830,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 3435, 3, 638, 866,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 3543, 3, 666, 902,
                                                                       ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 3651, 3, 794,
                                                                       1028, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 3786, 3, 830,
                                                                       1073, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 3921, 3, 866,
                                                                       1118, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 4056, 3, 1028,
                                                                       1273, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 4221, 3, 1073,
                                                                       1328, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 4386, 3, 1273,
                                                                       1515, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4584, 3, 8, 9,
                                                                       1581, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4590, 3, 9, 10,
                                                                       1584, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4596, 3, 10, 11,
                                                                       1587, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4602, 3, 11, 12,
                                                                       1590, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4608, 3, 12, 13,
                                                                       1593, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4614, 3, 13, 14,
                                                                       1596, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4620, 3, 14, 15,
                                                                       1599, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4626, 3, 15, 16,
                                                                       1602, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4632, 3, 16, 17,
                                                                       1605, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4638, 3, 17, 18,
                                                                       1608, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4644, 3, 18, 19,
                                                                       1611, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4650, 0, 3, 4584,
                                                                       1581, 4590, 1614, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4668, 0, 3, 4590,
                                                                       1584, 4596, 1623, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4686, 0, 3, 4596,
                                                                       1587, 4602, 1632, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4704, 0, 3, 4602,
                                                                       1590, 4608, 1641, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4722, 0, 3, 4608,
                                                                       1593, 4614, 1650, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4740, 0, 3, 4614,
                                                                       1596, 4620, 1659, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4758, 0, 3, 4620,
                                                                       1599, 4626, 1668, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4776, 0, 3, 4626,
                                                                       1602, 4632, 1677, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4794, 0, 3, 4632,
                                                                       1605, 4638, 1686, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4812, 0, 3, 4638,
                                                                       1608, 4644, 1695, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4830, 0, 3, 4650,
                                                                       1614, 4668, 57, 63, 1704,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4866, 0, 3, 4668,
                                                                       1623, 4686, 63, 69, 1722,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4902, 0, 3, 4686,
                                                                       1632, 4704, 69, 75, 1740,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4938, 0, 3, 4704,
                                                                       1641, 4722, 75, 81, 1758,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4974, 0, 3, 4722,
                                                                       1650, 4740, 81, 87, 1776,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5010, 0, 3, 4740,
                                                                       1659, 4758, 87, 93, 1794,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5046, 0, 3, 4758,
                                                                       1668, 4776, 93, 99, 1812,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5082, 0, 3, 4776,
                                                                       1677, 4794, 99, 105, 1830,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5118, 0, 3, 4794,
                                                                       1686, 4812, 105, 111,
                                                                       1848, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5154, 0, 3, 4830,
                                                                       1704, 4866, 123, 133,
                                                                       1866, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5214, 0, 3, 4866,
                                                                       1722, 4902, 133, 143,
                                                                       1896, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5274, 0, 3, 4902,
                                                                       1740, 4938, 143, 153,
                                                                       1926, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5334, 0, 3, 4938,
                                                                       1758, 4974, 153, 163,
                                                                       1956, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5394, 0, 3, 4974,
                                                                       1776, 5010, 163, 173,
                                                                       1986, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5454, 0, 3, 5010,
                                                                       1794, 5046, 173, 183,
                                                                       2016, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5514, 0, 3, 5046,
                                                                       1812, 5082, 183, 193,
                                                                       2046, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5574, 0, 3, 5082,
                                                                       1830, 5118, 193, 203,
                                                                       2076, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5634, 0, 3, 5154,
                                                                       1866, 5214, 223, 238,
                                                                       2106, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5724, 0, 3, 5214,
                                                                       1896, 5274, 238, 253,
                                                                       2151, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5814, 0, 3, 5274,
                                                                       1926, 5334, 253, 268,
                                                                       2196, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5904, 0, 3, 5334,
                                                                       1956, 5394, 268, 283,
                                                                       2241, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5994, 0, 3, 5394,
                                                                       1986, 5454, 283, 298,
                                                                       2286, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6084, 0, 3, 5454,
                                                                       2016, 5514, 298, 313,
                                                                       2331, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6174, 0, 3, 5514,
                                                                       2046, 5574, 313, 328,
                                                                       2376, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 6264, 0, 3, 5634,
                                                                       2106, 5724, 358, 379,
                                                                       2421, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 6390, 0, 3, 5724,
                                                                       2151, 5814, 379, 400,
                                                                       2484, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 6516, 0, 3, 5814,
                                                                       2196, 5904, 400, 421,
                                                                       2547, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 6642, 0, 3, 5904,
                                                                       2241, 5994, 421, 442,
                                                                       2610, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 6768, 0, 3, 5994,
                                                                       2286, 6084, 442, 463,
                                                                       2673, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 6894, 0, 3, 6084,
                                                                       2331, 6174, 463, 484,
                                                                       2736, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 7020, 0, 3, 6264,
                                                                       2421, 6390, 526, 554,
                                                                       2799, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 7188, 0, 3, 6390,
                                                                       2484, 6516, 554, 582,
                                                                       2883, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 7356, 0, 3, 6516,
                                                                       2547, 6642, 582, 610,
                                                                       2967, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 7524, 0, 3, 6642,
                                                                       2610, 6768, 610, 638,
                                                                       3051, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 7692, 0, 3, 6768,
                                                                       2673, 6894, 638, 666,
                                                                       3135, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 7860, 0, 3, 7020,
                                                                       2799, 7188, 722, 758,
                                                                       3219, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 8076, 0, 3, 7188,
                                                                       2883, 7356, 758, 794,
                                                                       3327, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 8292, 0, 3, 7356,
                                                                       2967, 7524, 794, 830,
                                                                       3435, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 8508, 0, 3, 7524,
                                                                       3051, 7692, 830, 866,
                                                                       3543, ncols, gamma, p,
                                                                       q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 8724, 0, 3, 7860,
                                                                       3219, 8076, 938, 983,
                                                                       3651, ncols, gamma, p,
                                                                       q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 8994, 0, 3, 8076,
                                                                       3327, 8292, 983, 1028,
                                                                       3786, ncols, gamma, p,
                                                                       q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 9264, 0, 3, 8292,
                                                                       3435, 8508, 1028, 1073,
                                                                       3921, ncols, gamma, p,
                                                                       q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 9534, 0, 3, 8724,
                                                                       3651, 8994, 1163, 1218,
                                                                       4056, ncols, gamma, p,
                                                                       q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 9864, 0, 3, 8994,
                                                                       3786, 9264, 1218, 1273,
                                                                       4221, ncols, gamma, p,
                                                                       q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 10194, 0, 3, 9534,
                                                                       4056, 9864, 1383, 1449,
                                                                       4386, ncols, gamma, p,
                                                                       q);

                    simdfunc::contract_primitives(buffer, 10590, 6264, 126, ncols);

                    simdfunc::contract_primitives(buffer, 10821, 7020, 168, ncols);

                    simdfunc::contract_primitives(buffer, 11129, 7860, 216, ncols);

                    simdfunc::contract_primitives(buffer, 11525, 8724, 270, ncols);

                    simdfunc::contract_primitives(buffer, 12020, 9534, 330, ncols);

                    simdfunc::contract_primitives(buffer, 12625, 10194, 396, ncols);
                }
            }
        }

        simdtrf::transform_d_inner(buffer, 10716, 10590, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 10989, 10821, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 11345, 11129, 36, 1, nmax);

        simdtrf::transform_d_inner(buffer, 11795, 11525, 45, 1, nmax);

        simdtrf::transform_d_inner(buffer, 12350, 12020, 55, 1, nmax);

        simdtrf::transform_d_inner(buffer, 13021, 12625, 66, 1, nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 13351, 10716, 10989, 5, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 13666, 10989, 11345, 5, nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 14086, 11345, 11795, 5, nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 14626, 11795, 12350, 5, nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 15301, 12350, 13021, 5, nmax);

        simdtrf::compute_hrr_hd_out_of_first(buffer, coordinates, 16126, 13351, 13666, 5, nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 16756, 13666, 14086, 5, nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 17596, 14086, 14626, 5, nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 18676, 14626, 15301, 5, nmax);

        simdtrf::compute_hrr_hf_out_of_first(buffer, coordinates, 20026, 16126, 16756, 5, nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 21076, 16756, 17596, 5, nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 22476, 17596, 18676, 5, nmax);

        simdtrf::compute_hrr_hg_out_of_first(buffer, coordinates, 24276, 20026, 21076, 5, nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 25851, 21076, 22476, 5, nmax);

        simdtrf::compute_hrr_hh_out_of_first(buffer, coordinates, 27951, 24276, 25851, 5, nmax);

        simdtrf::transform_h_inner(buffer, 30156, 27951, 21, 5, nmax);

        simdtrf::transform_h_outer(values + n * npairs, nvalues, buffer, 30156, 55, nmax);
    }

    for (size_t m = 0; m < 605; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
