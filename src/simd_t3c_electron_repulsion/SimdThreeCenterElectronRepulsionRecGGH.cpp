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


#include "SimdThreeCenterElectronRepulsionRecGGH.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferGD.hpp"
#include "SimdTransferGF.hpp"
#include "SimdTransferGG.hpp"
#include "SimdTransferGP.hpp"
#include "SimdTransferHD.hpp"
#include "SimdTransferHF.hpp"
#include "SimdTransferHP.hpp"
#include "SimdTransferID.hpp"
#include "SimdTransferIP.hpp"
#include "SimdTransferKP.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformH.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_ggh_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_ggh_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 55732, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 891 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 55732, 35648, 4145, dimensions);

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

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 7, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9, 10, 11, 12, 13}, ncols,
                                                        fj, 6, fq);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1163, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1166, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1169, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1172, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1175, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1178, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1181, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1184, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1187, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1190, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1193, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1196, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1199, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1202, 3, 8, 21,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1211, 3, 9, 24,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1220, 3, 10, 27,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1229, 3, 11, 30,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1238, 3, 12, 33,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1247, 3, 13, 36,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1256, 3, 14, 39,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1265, 3, 15, 42,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1274, 3, 16, 45,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1283, 3, 17, 48,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1292, 3, 18, 51,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1301, 3, 19, 54,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1310, 3, 21, 57,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1328, 3, 24, 63,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1346, 3, 27, 69,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1364, 3, 30, 75,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1382, 3, 33, 81,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1400, 3, 36, 87,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1418, 3, 39, 93,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1436, 3, 42, 99,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1454, 3, 45, 105,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1472, 3, 48, 111,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1490, 3, 51, 117,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1508, 3, 57, 123,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1538, 3, 63, 133,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1568, 3, 69, 143,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1598, 3, 75, 153,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1628, 3, 81, 163,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1658, 3, 87, 173,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1688, 3, 93, 183,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1718, 3, 99, 193,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1748, 3, 105, 203,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1778, 3, 111, 213,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1808, 3, 123, 223,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1853, 3, 133, 238,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1898, 3, 143, 253,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1943, 3, 153, 268,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1988, 3, 163, 283,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2033, 3, 173, 298,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2078, 3, 183, 313,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2123, 3, 193, 328,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2168, 3, 203, 343,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2213, 3, 223, 358,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2276, 3, 238, 379,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2339, 3, 253, 400,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2402, 3, 268, 421,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2465, 3, 283, 442,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2528, 3, 298, 463,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2591, 3, 313, 484,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2654, 3, 328, 505,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2717, 3, 358, 526,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2801, 3, 379, 554,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2885, 3, 400, 582,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2969, 3, 421, 610,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3053, 3, 442, 638,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3137, 3, 463, 666,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3221, 3, 484, 694,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 3305, 3, 526, 722,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 3413, 3, 554, 758,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 3521, 3, 582, 794,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 3629, 3, 610, 830,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 3737, 3, 638, 866,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 3845, 3, 666, 902,
                                                                       ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 3953, 3, 722, 938,
                                                                       ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 4088, 3, 758, 983,
                                                                       ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 4223, 3, 794,
                                                                       1028, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 4358, 3, 830,
                                                                       1073, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 4493, 3, 866,
                                                                       1118, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4628, 3, 8, 9,
                                                                       1169, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4634, 3, 9, 10,
                                                                       1172, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4640, 3, 10, 11,
                                                                       1175, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4646, 3, 11, 12,
                                                                       1178, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4652, 3, 12, 13,
                                                                       1181, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4658, 3, 13, 14,
                                                                       1184, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4664, 3, 14, 15,
                                                                       1187, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4670, 3, 15, 16,
                                                                       1190, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4676, 3, 16, 17,
                                                                       1193, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4682, 3, 17, 18,
                                                                       1196, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4688, 3, 18, 19,
                                                                       1199, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4694, 0, 3, 4628,
                                                                       1169, 4634, 1220, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4712, 0, 3, 4634,
                                                                       1172, 4640, 1229, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4730, 0, 3, 4640,
                                                                       1175, 4646, 1238, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4748, 0, 3, 4646,
                                                                       1178, 4652, 1247, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4766, 0, 3, 4652,
                                                                       1181, 4658, 1256, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4784, 0, 3, 4658,
                                                                       1184, 4664, 1265, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4802, 0, 3, 4664,
                                                                       1187, 4670, 1274, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4820, 0, 3, 4670,
                                                                       1190, 4676, 1283, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4838, 0, 3, 4676,
                                                                       1193, 4682, 1292, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4856, 0, 3, 4682,
                                                                       1196, 4688, 1301, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4874, 0, 3, 4694,
                                                                       1220, 4712, 57, 63, 1346,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4910, 0, 3, 4712,
                                                                       1229, 4730, 63, 69, 1364,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4946, 0, 3, 4730,
                                                                       1238, 4748, 69, 75, 1382,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4982, 0, 3, 4748,
                                                                       1247, 4766, 75, 81, 1400,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5018, 0, 3, 4766,
                                                                       1256, 4784, 81, 87, 1418,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5054, 0, 3, 4784,
                                                                       1265, 4802, 87, 93, 1436,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5090, 0, 3, 4802,
                                                                       1274, 4820, 93, 99, 1454,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5126, 0, 3, 4820,
                                                                       1283, 4838, 99, 105, 1472,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5162, 0, 3, 4838,
                                                                       1292, 4856, 105, 111,
                                                                       1490, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5198, 0, 3, 4874,
                                                                       1346, 4910, 123, 133,
                                                                       1568, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5258, 0, 3, 4910,
                                                                       1364, 4946, 133, 143,
                                                                       1598, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5318, 0, 3, 4946,
                                                                       1382, 4982, 143, 153,
                                                                       1628, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5378, 0, 3, 4982,
                                                                       1400, 5018, 153, 163,
                                                                       1658, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5438, 0, 3, 5018,
                                                                       1418, 5054, 163, 173,
                                                                       1688, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5498, 0, 3, 5054,
                                                                       1436, 5090, 173, 183,
                                                                       1718, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5558, 0, 3, 5090,
                                                                       1454, 5126, 183, 193,
                                                                       1748, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5618, 0, 3, 5126,
                                                                       1472, 5162, 193, 203,
                                                                       1778, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5678, 0, 3, 5198,
                                                                       1568, 5258, 223, 238,
                                                                       1898, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5768, 0, 3, 5258,
                                                                       1598, 5318, 238, 253,
                                                                       1943, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5858, 0, 3, 5318,
                                                                       1628, 5378, 253, 268,
                                                                       1988, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5948, 0, 3, 5378,
                                                                       1658, 5438, 268, 283,
                                                                       2033, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6038, 0, 3, 5438,
                                                                       1688, 5498, 283, 298,
                                                                       2078, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6128, 0, 3, 5498,
                                                                       1718, 5558, 298, 313,
                                                                       2123, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6218, 0, 3, 5558,
                                                                       1748, 5618, 313, 328,
                                                                       2168, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 6308, 0, 3, 5678,
                                                                       1898, 5768, 358, 379,
                                                                       2339, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 6434, 0, 3, 5768,
                                                                       1943, 5858, 379, 400,
                                                                       2402, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 6560, 0, 3, 5858,
                                                                       1988, 5948, 400, 421,
                                                                       2465, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 6686, 0, 3, 5948,
                                                                       2033, 6038, 421, 442,
                                                                       2528, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 6812, 0, 3, 6038,
                                                                       2078, 6128, 442, 463,
                                                                       2591, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 6938, 0, 3, 6128,
                                                                       2123, 6218, 463, 484,
                                                                       2654, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 7064, 0, 3, 6308,
                                                                       2339, 6434, 526, 554,
                                                                       2885, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 7232, 0, 3, 6434,
                                                                       2402, 6560, 554, 582,
                                                                       2969, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 7400, 0, 3, 6560,
                                                                       2465, 6686, 582, 610,
                                                                       3053, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 7568, 0, 3, 6686,
                                                                       2528, 6812, 610, 638,
                                                                       3137, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 7736, 0, 3, 6812,
                                                                       2591, 6938, 638, 666,
                                                                       3221, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 7904, 0, 3, 7064,
                                                                       2885, 7232, 722, 758,
                                                                       3521, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 8120, 0, 3, 7232,
                                                                       2969, 7400, 758, 794,
                                                                       3629, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 8336, 0, 3, 7400,
                                                                       3053, 7568, 794, 830,
                                                                       3737, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 8552, 0, 3, 7568,
                                                                       3137, 7736, 830, 866,
                                                                       3845, ncols, gamma, p,
                                                                       q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 8768, 0, 3, 7904,
                                                                       3521, 8120, 938, 983,
                                                                       4223, ncols, gamma, p,
                                                                       q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 9038, 0, 3, 8120,
                                                                       3629, 8336, 983, 1028,
                                                                       4358, ncols, gamma, p,
                                                                       q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 9308, 0, 3, 8336,
                                                                       3737, 8552, 1028, 1073,
                                                                       4493, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9578, 3, 1163,
                                                                       1166, 4628, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9588, 3, 1166,
                                                                       1169, 4634, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9598, 3, 1169,
                                                                       1172, 4640, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9608, 3, 1172,
                                                                       1175, 4646, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9618, 3, 1175,
                                                                       1178, 4652, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9628, 3, 1178,
                                                                       1181, 4658, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9638, 3, 1181,
                                                                       1184, 4664, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9648, 3, 1184,
                                                                       1187, 4670, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9658, 3, 1187,
                                                                       1190, 4676, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9668, 3, 1190,
                                                                       1193, 4682, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9678, 3, 1193,
                                                                       1196, 4688, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9688, 0, 3, 9578,
                                                                       4628, 9588, 1202, 1211,
                                                                       4694, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9718, 0, 3, 9588,
                                                                       4634, 9598, 1211, 1220,
                                                                       4712, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9748, 0, 3, 9598,
                                                                       4640, 9608, 1220, 1229,
                                                                       4730, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9778, 0, 3, 9608,
                                                                       4646, 9618, 1229, 1238,
                                                                       4748, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9808, 0, 3, 9618,
                                                                       4652, 9628, 1238, 1247,
                                                                       4766, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9838, 0, 3, 9628,
                                                                       4658, 9638, 1247, 1256,
                                                                       4784, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9868, 0, 3, 9638,
                                                                       4664, 9648, 1256, 1265,
                                                                       4802, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9898, 0, 3, 9648,
                                                                       4670, 9658, 1265, 1274,
                                                                       4820, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9928, 0, 3, 9658,
                                                                       4676, 9668, 1274, 1283,
                                                                       4838, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9958, 0, 3, 9668,
                                                                       4682, 9678, 1283, 1292,
                                                                       4856, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9988, 0, 3, 9688,
                                                                       4694, 9718, 1310, 1328,
                                                                       4874, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 10048, 0, 3, 9718,
                                                                       4712, 9748, 1328, 1346,
                                                                       4910, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 10108, 0, 3, 9748,
                                                                       4730, 9778, 1346, 1364,
                                                                       4946, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 10168, 0, 3, 9778,
                                                                       4748, 9808, 1364, 1382,
                                                                       4982, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 10228, 0, 3, 9808,
                                                                       4766, 9838, 1382, 1400,
                                                                       5018, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 10288, 0, 3, 9838,
                                                                       4784, 9868, 1400, 1418,
                                                                       5054, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 10348, 0, 3, 9868,
                                                                       4802, 9898, 1418, 1436,
                                                                       5090, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 10408, 0, 3, 9898,
                                                                       4820, 9928, 1436, 1454,
                                                                       5126, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 10468, 0, 3, 9928,
                                                                       4838, 9958, 1454, 1472,
                                                                       5162, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 10528, 0, 3, 9988,
                                                                       4874, 10048, 1508, 1538,
                                                                       5198, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 10628, 0, 3,
                                                                       10048, 4910, 10108, 1538,
                                                                       1568, 5258, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 10728, 0, 3,
                                                                       10108, 4946, 10168, 1568,
                                                                       1598, 5318, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 10828, 0, 3,
                                                                       10168, 4982, 10228, 1598,
                                                                       1628, 5378, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 10928, 0, 3,
                                                                       10228, 5018, 10288, 1628,
                                                                       1658, 5438, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 11028, 0, 3,
                                                                       10288, 5054, 10348, 1658,
                                                                       1688, 5498, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 11128, 0, 3,
                                                                       10348, 5090, 10408, 1688,
                                                                       1718, 5558, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 11228, 0, 3,
                                                                       10408, 5126, 10468, 1718,
                                                                       1748, 5618, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 11328, 0, 3,
                                                                       10528, 5198, 10628, 1808,
                                                                       1853, 5678, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 11478, 0, 3,
                                                                       10628, 5258, 10728, 1853,
                                                                       1898, 5768, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 11628, 0, 3,
                                                                       10728, 5318, 10828, 1898,
                                                                       1943, 5858, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 11778, 0, 3,
                                                                       10828, 5378, 10928, 1943,
                                                                       1988, 5948, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 11928, 0, 3,
                                                                       10928, 5438, 11028, 1988,
                                                                       2033, 6038, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 12078, 0, 3,
                                                                       11028, 5498, 11128, 2033,
                                                                       2078, 6128, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 12228, 0, 3,
                                                                       11128, 5558, 11228, 2078,
                                                                       2123, 6218, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 12378, 0, 3,
                                                                       11328, 5678, 11478, 2213,
                                                                       2276, 6308, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 12588, 0, 3,
                                                                       11478, 5768, 11628, 2276,
                                                                       2339, 6434, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 12798, 0, 3,
                                                                       11628, 5858, 11778, 2339,
                                                                       2402, 6560, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 13008, 0, 3,
                                                                       11778, 5948, 11928, 2402,
                                                                       2465, 6686, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 13218, 0, 3,
                                                                       11928, 6038, 12078, 2465,
                                                                       2528, 6812, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 13428, 0, 3,
                                                                       12078, 6128, 12228, 2528,
                                                                       2591, 6938, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 13638, 0, 3,
                                                                       12378, 6308, 12588, 2717,
                                                                       2801, 7064, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 13918, 0, 3,
                                                                       12588, 6434, 12798, 2801,
                                                                       2885, 7232, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 14198, 0, 3,
                                                                       12798, 6560, 13008, 2885,
                                                                       2969, 7400, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 14478, 0, 3,
                                                                       13008, 6686, 13218, 2969,
                                                                       3053, 7568, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 14758, 0, 3,
                                                                       13218, 6812, 13428, 3053,
                                                                       3137, 7736, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 15038, 0, 3,
                                                                       13638, 7064, 13918, 3305,
                                                                       3413, 7904, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 15398, 0, 3,
                                                                       13918, 7232, 14198, 3413,
                                                                       3521, 8120, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 15758, 0, 3,
                                                                       14198, 7400, 14478, 3521,
                                                                       3629, 8336, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 16118, 0, 3,
                                                                       14478, 7568, 14758, 3629,
                                                                       3737, 8552, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 16478, 0, 3,
                                                                       15038, 7904, 15398, 3953,
                                                                       4088, 8768, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 16928, 0, 3,
                                                                       15398, 8120, 15758, 4088,
                                                                       4223, 9038, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 17378, 0, 3,
                                                                       15758, 8336, 16118, 4223,
                                                                       4358, 9308, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 17828, 3, 4628,
                                                                       4634, 9598, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 17843, 3, 4634,
                                                                       4640, 9608, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 17858, 3, 4640,
                                                                       4646, 9618, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 17873, 3, 4646,
                                                                       4652, 9628, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 17888, 3, 4652,
                                                                       4658, 9638, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 17903, 3, 4658,
                                                                       4664, 9648, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 17918, 3, 4664,
                                                                       4670, 9658, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 17933, 3, 4670,
                                                                       4676, 9668, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 17948, 3, 4676,
                                                                       4682, 9678, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 17963, 0, 3,
                                                                       17828, 9598, 17843, 4694,
                                                                       4712, 9748, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 18008, 0, 3,
                                                                       17843, 9608, 17858, 4712,
                                                                       4730, 9778, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 18053, 0, 3,
                                                                       17858, 9618, 17873, 4730,
                                                                       4748, 9808, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 18098, 0, 3,
                                                                       17873, 9628, 17888, 4748,
                                                                       4766, 9838, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 18143, 0, 3,
                                                                       17888, 9638, 17903, 4766,
                                                                       4784, 9868, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 18188, 0, 3,
                                                                       17903, 9648, 17918, 4784,
                                                                       4802, 9898, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 18233, 0, 3,
                                                                       17918, 9658, 17933, 4802,
                                                                       4820, 9928, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 18278, 0, 3,
                                                                       17933, 9668, 17948, 4820,
                                                                       4838, 9958, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 18323, 0, 3,
                                                                       17963, 9748, 18008, 4874,
                                                                       4910, 10108, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 18413, 0, 3,
                                                                       18008, 9778, 18053, 4910,
                                                                       4946, 10168, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 18503, 0, 3,
                                                                       18053, 9808, 18098, 4946,
                                                                       4982, 10228, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 18593, 0, 3,
                                                                       18098, 9838, 18143, 4982,
                                                                       5018, 10288, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 18683, 0, 3,
                                                                       18143, 9868, 18188, 5018,
                                                                       5054, 10348, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 18773, 0, 3,
                                                                       18188, 9898, 18233, 5054,
                                                                       5090, 10408, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 18863, 0, 3,
                                                                       18233, 9928, 18278, 5090,
                                                                       5126, 10468, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 18953, 0, 3,
                                                                       18323, 10108, 18413, 5198,
                                                                       5258, 10728, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 19103, 0, 3,
                                                                       18413, 10168, 18503, 5258,
                                                                       5318, 10828, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 19253, 0, 3,
                                                                       18503, 10228, 18593, 5318,
                                                                       5378, 10928, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 19403, 0, 3,
                                                                       18593, 10288, 18683, 5378,
                                                                       5438, 11028, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 19553, 0, 3,
                                                                       18683, 10348, 18773, 5438,
                                                                       5498, 11128, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 19703, 0, 3,
                                                                       18773, 10408, 18863, 5498,
                                                                       5558, 11228, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 19853, 0, 3,
                                                                       18953, 10728, 19103, 5678,
                                                                       5768, 11628, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 20078, 0, 3,
                                                                       19103, 10828, 19253, 5768,
                                                                       5858, 11778, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 20303, 0, 3,
                                                                       19253, 10928, 19403, 5858,
                                                                       5948, 11928, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 20528, 0, 3,
                                                                       19403, 11028, 19553, 5948,
                                                                       6038, 12078, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 20753, 0, 3,
                                                                       19553, 11128, 19703, 6038,
                                                                       6128, 12228, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 20978, 0, 3,
                                                                       19853, 11628, 20078, 6308,
                                                                       6434, 12798, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 21293, 0, 3,
                                                                       20078, 11778, 20303, 6434,
                                                                       6560, 13008, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 21608, 0, 3,
                                                                       20303, 11928, 20528, 6560,
                                                                       6686, 13218, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 21923, 0, 3,
                                                                       20528, 12078, 20753, 6686,
                                                                       6812, 13428, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 22238, 0, 3,
                                                                       20978, 12798, 21293, 7064,
                                                                       7232, 14198, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 22658, 0, 3,
                                                                       21293, 13008, 21608, 7232,
                                                                       7400, 14478, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 23078, 0, 3,
                                                                       21608, 13218, 21923, 7400,
                                                                       7568, 14758, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 23498, 0, 3,
                                                                       22238, 14198, 22658, 7904,
                                                                       8120, 15758, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 24038, 0, 3,
                                                                       22658, 14478, 23078, 8120,
                                                                       8336, 16118, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 24578, 0, 3,
                                                                       23498, 15758, 24038, 8768,
                                                                       9038, 17378, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25253, 3, 9578,
                                                                       9588, 17828, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25274, 3, 9588,
                                                                       9598, 17843, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25295, 3, 9598,
                                                                       9608, 17858, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25316, 3, 9608,
                                                                       9618, 17873, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25337, 3, 9618,
                                                                       9628, 17888, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25358, 3, 9628,
                                                                       9638, 17903, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25379, 3, 9638,
                                                                       9648, 17918, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25400, 3, 9648,
                                                                       9658, 17933, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25421, 3, 9658,
                                                                       9668, 17948, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 25442, 0, 3,
                                                                       25253, 17828, 25274, 9688,
                                                                       9718, 17963, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 25505, 0, 3,
                                                                       25274, 17843, 25295, 9718,
                                                                       9748, 18008, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 25568, 0, 3,
                                                                       25295, 17858, 25316, 9748,
                                                                       9778, 18053, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 25631, 0, 3,
                                                                       25316, 17873, 25337, 9778,
                                                                       9808, 18098, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 25694, 0, 3,
                                                                       25337, 17888, 25358, 9808,
                                                                       9838, 18143, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 25757, 0, 3,
                                                                       25358, 17903, 25379, 9838,
                                                                       9868, 18188, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 25820, 0, 3,
                                                                       25379, 17918, 25400, 9868,
                                                                       9898, 18233, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 25883, 0, 3,
                                                                       25400, 17933, 25421, 9898,
                                                                       9928, 18278, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 25946, 0, 3,
                                                                       25442, 17963, 25505, 9988,
                                                                       10048, 18323, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 26072, 0, 3,
                                                                       25505, 18008, 25568,
                                                                       10048, 10108, 18413,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 26198, 0, 3,
                                                                       25568, 18053, 25631,
                                                                       10108, 10168, 18503,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 26324, 0, 3,
                                                                       25631, 18098, 25694,
                                                                       10168, 10228, 18593,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 26450, 0, 3,
                                                                       25694, 18143, 25757,
                                                                       10228, 10288, 18683,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 26576, 0, 3,
                                                                       25757, 18188, 25820,
                                                                       10288, 10348, 18773,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 26702, 0, 3,
                                                                       25820, 18233, 25883,
                                                                       10348, 10408, 18863,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 26828, 0, 3,
                                                                       25946, 18323, 26072,
                                                                       10528, 10628, 18953,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 27038, 0, 3,
                                                                       26072, 18413, 26198,
                                                                       10628, 10728, 19103,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 27248, 0, 3,
                                                                       26198, 18503, 26324,
                                                                       10728, 10828, 19253,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 27458, 0, 3,
                                                                       26324, 18593, 26450,
                                                                       10828, 10928, 19403,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 27668, 0, 3,
                                                                       26450, 18683, 26576,
                                                                       10928, 11028, 19553,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 27878, 0, 3,
                                                                       26576, 18773, 26702,
                                                                       11028, 11128, 19703,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 28088, 0, 3,
                                                                       26828, 18953, 27038,
                                                                       11328, 11478, 19853,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 28403, 0, 3,
                                                                       27038, 19103, 27248,
                                                                       11478, 11628, 20078,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 28718, 0, 3,
                                                                       27248, 19253, 27458,
                                                                       11628, 11778, 20303,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 29033, 0, 3,
                                                                       27458, 19403, 27668,
                                                                       11778, 11928, 20528,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 29348, 0, 3,
                                                                       27668, 19553, 27878,
                                                                       11928, 12078, 20753,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 29663, 0, 3,
                                                                       28088, 19853, 28403,
                                                                       12378, 12588, 20978,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 30104, 0, 3,
                                                                       28403, 20078, 28718,
                                                                       12588, 12798, 21293,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 30545, 0, 3,
                                                                       28718, 20303, 29033,
                                                                       12798, 13008, 21608,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 30986, 0, 3,
                                                                       29033, 20528, 29348,
                                                                       13008, 13218, 21923,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 31427, 0, 3,
                                                                       29663, 20978, 30104,
                                                                       13638, 13918, 22238,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 32015, 0, 3,
                                                                       30104, 21293, 30545,
                                                                       13918, 14198, 22658,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 32603, 0, 3,
                                                                       30545, 21608, 30986,
                                                                       14198, 14478, 23078,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 33191, 0, 3,
                                                                       31427, 22238, 32015,
                                                                       15038, 15398, 23498,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 33947, 0, 3,
                                                                       32015, 22658, 32603,
                                                                       15398, 15758, 24038,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 34703, 0, 3,
                                                                       33191, 23498, 33947,
                                                                       16478, 16928, 24578,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 35648, 28088, 315, ncols);

                    simdfunc::contract_primitives(buffer, 36128, 29663, 441, ncols);

                    simdfunc::contract_primitives(buffer, 36800, 31427, 588, ncols);

                    simdfunc::contract_primitives(buffer, 37696, 33191, 756, ncols);

                    simdfunc::contract_primitives(buffer, 38848, 34703, 945, ncols);
                }
            }
        }

        simdtrf::transform_h_inner(buffer, 35963, 35648, 15, 1, nmax);

        simdtrf::transform_h_inner(buffer, 36569, 36128, 21, 1, nmax);

        simdtrf::transform_h_inner(buffer, 37388, 36800, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 38452, 37696, 36, 1, nmax);

        simdtrf::transform_h_inner(buffer, 39793, 38848, 45, 1, nmax);

        simdtrf::compute_hrr_gp_out_of_first(buffer, coordinates, 40288, 35963, 36569, 11,
                                             nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 40783, 36569, 37388, 11,
                                             nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 41476, 37388, 38452, 11,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 42400, 38452, 39793, 11,
                                             nmax);

        simdtrf::compute_hrr_gd_out_of_first(buffer, coordinates, 43588, 40288, 40783, 11,
                                             nmax);

        simdtrf::compute_hrr_hd_out_of_first(buffer, coordinates, 44578, 40783, 41476, 11,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 45964, 41476, 42400, 11,
                                             nmax);

        simdtrf::compute_hrr_gf_out_of_first(buffer, coordinates, 47812, 43588, 44578, 11,
                                             nmax);

        simdtrf::compute_hrr_hf_out_of_first(buffer, coordinates, 49462, 44578, 45964, 11,
                                             nmax);

        simdtrf::compute_hrr_gg_out_of_first(buffer, coordinates, 51772, 47812, 49462, 11,
                                             nmax);

        simdtrf::transform_g_inner(buffer, 54247, 51772, 15, 11, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 54247, 99, nmax);
    }

    for (size_t m = 0; m < 891; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
