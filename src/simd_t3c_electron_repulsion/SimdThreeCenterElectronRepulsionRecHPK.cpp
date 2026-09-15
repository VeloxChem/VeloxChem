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


#include "SimdThreeCenterElectronRepulsionRecHPK.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferHP.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformK.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_hpk_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_hpk_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 43037, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 495 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 43037, 38648, 2079, dimensions);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 722, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 725, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 728, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 731, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 734, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 737, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 740, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 743, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 746, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 749, 3, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 752, 3, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 755, 3, 19, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 758, 3, 20, ncols,
                                                                       p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 761, 3, 8, 21,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 770, 3, 9, 24,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 779, 3, 10, 27,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 788, 3, 11, 30,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 797, 3, 12, 33,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 806, 3, 13, 36,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 815, 3, 14, 39,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 824, 3, 15, 42,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 833, 3, 16, 45,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 842, 3, 17, 48,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 851, 3, 18, 51,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 860, 3, 19, 54,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 869, 3, 21, 57,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 887, 3, 24, 63,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 905, 3, 27, 69,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 923, 3, 30, 75,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 941, 3, 33, 81,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 959, 3, 36, 87,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 977, 3, 39, 93,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 995, 3, 42, 99,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1013, 3, 45, 105,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1031, 3, 48, 111,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1049, 3, 51, 117,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1067, 3, 57, 123,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1097, 3, 63, 133,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1127, 3, 69, 143,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1157, 3, 75, 153,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1187, 3, 81, 163,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1217, 3, 87, 173,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1247, 3, 93, 183,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1277, 3, 99, 193,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1307, 3, 105, 203,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1337, 3, 111, 213,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1367, 3, 123, 223,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1412, 3, 133, 238,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1457, 3, 143, 253,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1502, 3, 153, 268,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1547, 3, 163, 283,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1592, 3, 173, 298,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1637, 3, 183, 313,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1682, 3, 193, 328,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1727, 3, 203, 343,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1772, 3, 223, 358,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1835, 3, 238, 379,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1898, 3, 253, 400,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1961, 3, 268, 421,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2024, 3, 283, 442,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2087, 3, 298, 463,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2150, 3, 313, 484,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2213, 3, 328, 505,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2276, 3, 358, 526,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2360, 3, 379, 554,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2444, 3, 400, 582,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2528, 3, 421, 610,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2612, 3, 442, 638,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2696, 3, 463, 666,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2780, 3, 484, 694,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2864, 3, 8, 9,
                                                                       728, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2870, 3, 9, 10,
                                                                       731, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2876, 3, 10, 11,
                                                                       734, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2882, 3, 11, 12,
                                                                       737, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2888, 3, 12, 13,
                                                                       740, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2894, 3, 13, 14,
                                                                       743, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2900, 3, 14, 15,
                                                                       746, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2906, 3, 15, 16,
                                                                       749, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2912, 3, 16, 17,
                                                                       752, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2918, 3, 17, 18,
                                                                       755, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2924, 3, 18, 19,
                                                                       758, ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2930, 0, 3, 2864,
                                                                       728, 2870, 779, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2948, 0, 3, 2870,
                                                                       731, 2876, 788, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2966, 0, 3, 2876,
                                                                       734, 2882, 797, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2984, 0, 3, 2882,
                                                                       737, 2888, 806, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3002, 0, 3, 2888,
                                                                       740, 2894, 815, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3020, 0, 3, 2894,
                                                                       743, 2900, 824, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3038, 0, 3, 2900,
                                                                       746, 2906, 833, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3056, 0, 3, 2906,
                                                                       749, 2912, 842, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3074, 0, 3, 2912,
                                                                       752, 2918, 851, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3092, 0, 3, 2918,
                                                                       755, 2924, 860, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3110, 0, 3, 2930,
                                                                       779, 2948, 57, 63, 905,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3146, 0, 3, 2948,
                                                                       788, 2966, 63, 69, 923,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3182, 0, 3, 2966,
                                                                       797, 2984, 69, 75, 941,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3218, 0, 3, 2984,
                                                                       806, 3002, 75, 81, 959,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3254, 0, 3, 3002,
                                                                       815, 3020, 81, 87, 977,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3290, 0, 3, 3020,
                                                                       824, 3038, 87, 93, 995,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3326, 0, 3, 3038,
                                                                       833, 3056, 93, 99, 1013,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3362, 0, 3, 3056,
                                                                       842, 3074, 99, 105, 1031,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3398, 0, 3, 3074,
                                                                       851, 3092, 105, 111, 1049,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3434, 0, 3, 3110,
                                                                       905, 3146, 123, 133, 1127,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3494, 0, 3, 3146,
                                                                       923, 3182, 133, 143, 1157,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3554, 0, 3, 3182,
                                                                       941, 3218, 143, 153, 1187,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3614, 0, 3, 3218,
                                                                       959, 3254, 153, 163, 1217,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3674, 0, 3, 3254,
                                                                       977, 3290, 163, 173, 1247,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3734, 0, 3, 3290,
                                                                       995, 3326, 173, 183, 1277,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3794, 0, 3, 3326,
                                                                       1013, 3362, 183, 193,
                                                                       1307, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3854, 0, 3, 3362,
                                                                       1031, 3398, 193, 203,
                                                                       1337, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 3914, 0, 3, 3434,
                                                                       1127, 3494, 223, 238,
                                                                       1457, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4004, 0, 3, 3494,
                                                                       1157, 3554, 238, 253,
                                                                       1502, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4094, 0, 3, 3554,
                                                                       1187, 3614, 253, 268,
                                                                       1547, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4184, 0, 3, 3614,
                                                                       1217, 3674, 268, 283,
                                                                       1592, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4274, 0, 3, 3674,
                                                                       1247, 3734, 283, 298,
                                                                       1637, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4364, 0, 3, 3734,
                                                                       1277, 3794, 298, 313,
                                                                       1682, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4454, 0, 3, 3794,
                                                                       1307, 3854, 313, 328,
                                                                       1727, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 4544, 0, 3, 3914,
                                                                       1457, 4004, 358, 379,
                                                                       1898, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 4670, 0, 3, 4004,
                                                                       1502, 4094, 379, 400,
                                                                       1961, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 4796, 0, 3, 4094,
                                                                       1547, 4184, 400, 421,
                                                                       2024, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 4922, 0, 3, 4184,
                                                                       1592, 4274, 421, 442,
                                                                       2087, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 5048, 0, 3, 4274,
                                                                       1637, 4364, 442, 463,
                                                                       2150, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 5174, 0, 3, 4364,
                                                                       1682, 4454, 463, 484,
                                                                       2213, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 5300, 0, 3, 4544,
                                                                       1898, 4670, 526, 554,
                                                                       2444, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 5468, 0, 3, 4670,
                                                                       1961, 4796, 554, 582,
                                                                       2528, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 5636, 0, 3, 4796,
                                                                       2024, 4922, 582, 610,
                                                                       2612, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 5804, 0, 3, 4922,
                                                                       2087, 5048, 610, 638,
                                                                       2696, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 5972, 0, 3, 5048,
                                                                       2150, 5174, 638, 666,
                                                                       2780, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6140, 3, 722, 725,
                                                                       2864, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6150, 3, 725, 728,
                                                                       2870, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6160, 3, 728, 731,
                                                                       2876, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6170, 3, 731, 734,
                                                                       2882, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6180, 3, 734, 737,
                                                                       2888, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6190, 3, 737, 740,
                                                                       2894, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6200, 3, 740, 743,
                                                                       2900, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6210, 3, 743, 746,
                                                                       2906, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6220, 3, 746, 749,
                                                                       2912, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6230, 3, 749, 752,
                                                                       2918, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6240, 3, 752, 755,
                                                                       2924, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6250, 0, 3, 6140,
                                                                       2864, 6150, 761, 770,
                                                                       2930, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6280, 0, 3, 6150,
                                                                       2870, 6160, 770, 779,
                                                                       2948, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6310, 0, 3, 6160,
                                                                       2876, 6170, 779, 788,
                                                                       2966, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6340, 0, 3, 6170,
                                                                       2882, 6180, 788, 797,
                                                                       2984, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6370, 0, 3, 6180,
                                                                       2888, 6190, 797, 806,
                                                                       3002, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6400, 0, 3, 6190,
                                                                       2894, 6200, 806, 815,
                                                                       3020, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6430, 0, 3, 6200,
                                                                       2900, 6210, 815, 824,
                                                                       3038, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6460, 0, 3, 6210,
                                                                       2906, 6220, 824, 833,
                                                                       3056, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6490, 0, 3, 6220,
                                                                       2912, 6230, 833, 842,
                                                                       3074, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6520, 0, 3, 6230,
                                                                       2918, 6240, 842, 851,
                                                                       3092, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6550, 0, 3, 6250,
                                                                       2930, 6280, 869, 887,
                                                                       3110, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6610, 0, 3, 6280,
                                                                       2948, 6310, 887, 905,
                                                                       3146, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6670, 0, 3, 6310,
                                                                       2966, 6340, 905, 923,
                                                                       3182, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6730, 0, 3, 6340,
                                                                       2984, 6370, 923, 941,
                                                                       3218, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6790, 0, 3, 6370,
                                                                       3002, 6400, 941, 959,
                                                                       3254, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6850, 0, 3, 6400,
                                                                       3020, 6430, 959, 977,
                                                                       3290, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6910, 0, 3, 6430,
                                                                       3038, 6460, 977, 995,
                                                                       3326, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6970, 0, 3, 6460,
                                                                       3056, 6490, 995, 1013,
                                                                       3362, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 7030, 0, 3, 6490,
                                                                       3074, 6520, 1013, 1031,
                                                                       3398, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7090, 0, 3, 6550,
                                                                       3110, 6610, 1067, 1097,
                                                                       3434, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7190, 0, 3, 6610,
                                                                       3146, 6670, 1097, 1127,
                                                                       3494, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7290, 0, 3, 6670,
                                                                       3182, 6730, 1127, 1157,
                                                                       3554, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7390, 0, 3, 6730,
                                                                       3218, 6790, 1157, 1187,
                                                                       3614, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7490, 0, 3, 6790,
                                                                       3254, 6850, 1187, 1217,
                                                                       3674, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7590, 0, 3, 6850,
                                                                       3290, 6910, 1217, 1247,
                                                                       3734, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7690, 0, 3, 6910,
                                                                       3326, 6970, 1247, 1277,
                                                                       3794, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7790, 0, 3, 6970,
                                                                       3362, 7030, 1277, 1307,
                                                                       3854, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 7890, 0, 3, 7090,
                                                                       3434, 7190, 1367, 1412,
                                                                       3914, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 8040, 0, 3, 7190,
                                                                       3494, 7290, 1412, 1457,
                                                                       4004, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 8190, 0, 3, 7290,
                                                                       3554, 7390, 1457, 1502,
                                                                       4094, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 8340, 0, 3, 7390,
                                                                       3614, 7490, 1502, 1547,
                                                                       4184, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 8490, 0, 3, 7490,
                                                                       3674, 7590, 1547, 1592,
                                                                       4274, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 8640, 0, 3, 7590,
                                                                       3734, 7690, 1592, 1637,
                                                                       4364, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 8790, 0, 3, 7690,
                                                                       3794, 7790, 1637, 1682,
                                                                       4454, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 8940, 0, 3, 7890,
                                                                       3914, 8040, 1772, 1835,
                                                                       4544, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 9150, 0, 3, 8040,
                                                                       4004, 8190, 1835, 1898,
                                                                       4670, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 9360, 0, 3, 8190,
                                                                       4094, 8340, 1898, 1961,
                                                                       4796, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 9570, 0, 3, 8340,
                                                                       4184, 8490, 1961, 2024,
                                                                       4922, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 9780, 0, 3, 8490,
                                                                       4274, 8640, 2024, 2087,
                                                                       5048, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 9990, 0, 3, 8640,
                                                                       4364, 8790, 2087, 2150,
                                                                       5174, ncols, gamma, p,
                                                                       q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 10200, 0, 3, 8940,
                                                                       4544, 9150, 2276, 2360,
                                                                       5300, ncols, gamma, p,
                                                                       q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 10480, 0, 3, 9150,
                                                                       4670, 9360, 2360, 2444,
                                                                       5468, ncols, gamma, p,
                                                                       q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 10760, 0, 3, 9360,
                                                                       4796, 9570, 2444, 2528,
                                                                       5636, ncols, gamma, p,
                                                                       q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 11040, 0, 3, 9570,
                                                                       4922, 9780, 2528, 2612,
                                                                       5804, ncols, gamma, p,
                                                                       q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 11320, 0, 3, 9780,
                                                                       5048, 9990, 2612, 2696,
                                                                       5972, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 11600, 3, 2864,
                                                                       2870, 6160, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 11615, 3, 2870,
                                                                       2876, 6170, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 11630, 3, 2876,
                                                                       2882, 6180, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 11645, 3, 2882,
                                                                       2888, 6190, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 11660, 3, 2888,
                                                                       2894, 6200, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 11675, 3, 2894,
                                                                       2900, 6210, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 11690, 3, 2900,
                                                                       2906, 6220, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 11705, 3, 2906,
                                                                       2912, 6230, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 11720, 3, 2912,
                                                                       2918, 6240, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 11735, 0, 3,
                                                                       11600, 6160, 11615, 2930,
                                                                       2948, 6310, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 11780, 0, 3,
                                                                       11615, 6170, 11630, 2948,
                                                                       2966, 6340, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 11825, 0, 3,
                                                                       11630, 6180, 11645, 2966,
                                                                       2984, 6370, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 11870, 0, 3,
                                                                       11645, 6190, 11660, 2984,
                                                                       3002, 6400, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 11915, 0, 3,
                                                                       11660, 6200, 11675, 3002,
                                                                       3020, 6430, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 11960, 0, 3,
                                                                       11675, 6210, 11690, 3020,
                                                                       3038, 6460, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 12005, 0, 3,
                                                                       11690, 6220, 11705, 3038,
                                                                       3056, 6490, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 12050, 0, 3,
                                                                       11705, 6230, 11720, 3056,
                                                                       3074, 6520, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 12095, 0, 3,
                                                                       11735, 6310, 11780, 3110,
                                                                       3146, 6670, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 12185, 0, 3,
                                                                       11780, 6340, 11825, 3146,
                                                                       3182, 6730, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 12275, 0, 3,
                                                                       11825, 6370, 11870, 3182,
                                                                       3218, 6790, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 12365, 0, 3,
                                                                       11870, 6400, 11915, 3218,
                                                                       3254, 6850, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 12455, 0, 3,
                                                                       11915, 6430, 11960, 3254,
                                                                       3290, 6910, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 12545, 0, 3,
                                                                       11960, 6460, 12005, 3290,
                                                                       3326, 6970, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 12635, 0, 3,
                                                                       12005, 6490, 12050, 3326,
                                                                       3362, 7030, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 12725, 0, 3,
                                                                       12095, 6670, 12185, 3434,
                                                                       3494, 7290, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 12875, 0, 3,
                                                                       12185, 6730, 12275, 3494,
                                                                       3554, 7390, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 13025, 0, 3,
                                                                       12275, 6790, 12365, 3554,
                                                                       3614, 7490, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 13175, 0, 3,
                                                                       12365, 6850, 12455, 3614,
                                                                       3674, 7590, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 13325, 0, 3,
                                                                       12455, 6910, 12545, 3674,
                                                                       3734, 7690, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 13475, 0, 3,
                                                                       12545, 6970, 12635, 3734,
                                                                       3794, 7790, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 13625, 0, 3,
                                                                       12725, 7290, 12875, 3914,
                                                                       4004, 8190, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 13850, 0, 3,
                                                                       12875, 7390, 13025, 4004,
                                                                       4094, 8340, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 14075, 0, 3,
                                                                       13025, 7490, 13175, 4094,
                                                                       4184, 8490, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 14300, 0, 3,
                                                                       13175, 7590, 13325, 4184,
                                                                       4274, 8640, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 14525, 0, 3,
                                                                       13325, 7690, 13475, 4274,
                                                                       4364, 8790, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 14750, 0, 3,
                                                                       13625, 8190, 13850, 4544,
                                                                       4670, 9360, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 15065, 0, 3,
                                                                       13850, 8340, 14075, 4670,
                                                                       4796, 9570, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 15380, 0, 3,
                                                                       14075, 8490, 14300, 4796,
                                                                       4922, 9780, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 15695, 0, 3,
                                                                       14300, 8640, 14525, 4922,
                                                                       5048, 9990, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 16010, 0, 3,
                                                                       14750, 9360, 15065, 5300,
                                                                       5468, 10760, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 16430, 0, 3,
                                                                       15065, 9570, 15380, 5468,
                                                                       5636, 11040, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 16850, 0, 3,
                                                                       15380, 9780, 15695, 5636,
                                                                       5804, 11320, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17270, 3, 6140,
                                                                       6150, 11600, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17291, 3, 6150,
                                                                       6160, 11615, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17312, 3, 6160,
                                                                       6170, 11630, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17333, 3, 6170,
                                                                       6180, 11645, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17354, 3, 6180,
                                                                       6190, 11660, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17375, 3, 6190,
                                                                       6200, 11675, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17396, 3, 6200,
                                                                       6210, 11690, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17417, 3, 6210,
                                                                       6220, 11705, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17438, 3, 6220,
                                                                       6230, 11720, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 17459, 0, 3,
                                                                       17270, 11600, 17291, 6250,
                                                                       6280, 11735, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 17522, 0, 3,
                                                                       17291, 11615, 17312, 6280,
                                                                       6310, 11780, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 17585, 0, 3,
                                                                       17312, 11630, 17333, 6310,
                                                                       6340, 11825, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 17648, 0, 3,
                                                                       17333, 11645, 17354, 6340,
                                                                       6370, 11870, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 17711, 0, 3,
                                                                       17354, 11660, 17375, 6370,
                                                                       6400, 11915, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 17774, 0, 3,
                                                                       17375, 11675, 17396, 6400,
                                                                       6430, 11960, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 17837, 0, 3,
                                                                       17396, 11690, 17417, 6430,
                                                                       6460, 12005, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 17900, 0, 3,
                                                                       17417, 11705, 17438, 6460,
                                                                       6490, 12050, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 17963, 0, 3,
                                                                       17459, 11735, 17522, 6550,
                                                                       6610, 12095, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 18089, 0, 3,
                                                                       17522, 11780, 17585, 6610,
                                                                       6670, 12185, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 18215, 0, 3,
                                                                       17585, 11825, 17648, 6670,
                                                                       6730, 12275, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 18341, 0, 3,
                                                                       17648, 11870, 17711, 6730,
                                                                       6790, 12365, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 18467, 0, 3,
                                                                       17711, 11915, 17774, 6790,
                                                                       6850, 12455, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 18593, 0, 3,
                                                                       17774, 11960, 17837, 6850,
                                                                       6910, 12545, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 18719, 0, 3,
                                                                       17837, 12005, 17900, 6910,
                                                                       6970, 12635, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 18845, 0, 3,
                                                                       17963, 12095, 18089, 7090,
                                                                       7190, 12725, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 19055, 0, 3,
                                                                       18089, 12185, 18215, 7190,
                                                                       7290, 12875, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 19265, 0, 3,
                                                                       18215, 12275, 18341, 7290,
                                                                       7390, 13025, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 19475, 0, 3,
                                                                       18341, 12365, 18467, 7390,
                                                                       7490, 13175, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 19685, 0, 3,
                                                                       18467, 12455, 18593, 7490,
                                                                       7590, 13325, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 19895, 0, 3,
                                                                       18593, 12545, 18719, 7590,
                                                                       7690, 13475, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 20105, 0, 3,
                                                                       18845, 12725, 19055, 7890,
                                                                       8040, 13625, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 20420, 0, 3,
                                                                       19055, 12875, 19265, 8040,
                                                                       8190, 13850, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 20735, 0, 3,
                                                                       19265, 13025, 19475, 8190,
                                                                       8340, 14075, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 21050, 0, 3,
                                                                       19475, 13175, 19685, 8340,
                                                                       8490, 14300, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 21365, 0, 3,
                                                                       19685, 13325, 19895, 8490,
                                                                       8640, 14525, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 21680, 0, 3,
                                                                       20105, 13625, 20420, 8940,
                                                                       9150, 14750, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 22121, 0, 3,
                                                                       20420, 13850, 20735, 9150,
                                                                       9360, 15065, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 22562, 0, 3,
                                                                       20735, 14075, 21050, 9360,
                                                                       9570, 15380, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 23003, 0, 3,
                                                                       21050, 14300, 21365, 9570,
                                                                       9780, 15695, ncols, gamma,
                                                                       p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 23444, 0, 3,
                                                                       21680, 14750, 22121,
                                                                       10200, 10480, 16010,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 24032, 0, 3,
                                                                       22121, 15065, 22562,
                                                                       10480, 10760, 16430,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 24620, 0, 3,
                                                                       22562, 15380, 23003,
                                                                       10760, 11040, 16850,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 25208, 3, 11600,
                                                                       11615, 17312, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 25236, 3, 11615,
                                                                       11630, 17333, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 25264, 3, 11630,
                                                                       11645, 17354, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 25292, 3, 11645,
                                                                       11660, 17375, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 25320, 3, 11660,
                                                                       11675, 17396, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 25348, 3, 11675,
                                                                       11690, 17417, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 25376, 3, 11690,
                                                                       11705, 17438, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 25404, 0, 3,
                                                                       25208, 17312, 25236,
                                                                       11735, 11780, 17585,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 25488, 0, 3,
                                                                       25236, 17333, 25264,
                                                                       11780, 11825, 17648,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 25572, 0, 3,
                                                                       25264, 17354, 25292,
                                                                       11825, 11870, 17711,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 25656, 0, 3,
                                                                       25292, 17375, 25320,
                                                                       11870, 11915, 17774,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 25740, 0, 3,
                                                                       25320, 17396, 25348,
                                                                       11915, 11960, 17837,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 25824, 0, 3,
                                                                       25348, 17417, 25376,
                                                                       11960, 12005, 17900,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 25908, 0, 3,
                                                                       25404, 17585, 25488,
                                                                       12095, 12185, 18215,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 26076, 0, 3,
                                                                       25488, 17648, 25572,
                                                                       12185, 12275, 18341,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 26244, 0, 3,
                                                                       25572, 17711, 25656,
                                                                       12275, 12365, 18467,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 26412, 0, 3,
                                                                       25656, 17774, 25740,
                                                                       12365, 12455, 18593,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 26580, 0, 3,
                                                                       25740, 17837, 25824,
                                                                       12455, 12545, 18719,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 26748, 0, 3,
                                                                       25908, 18215, 26076,
                                                                       12725, 12875, 19265,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 27028, 0, 3,
                                                                       26076, 18341, 26244,
                                                                       12875, 13025, 19475,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 27308, 0, 3,
                                                                       26244, 18467, 26412,
                                                                       13025, 13175, 19685,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 27588, 0, 3,
                                                                       26412, 18593, 26580,
                                                                       13175, 13325, 19895,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 27868, 0, 3,
                                                                       26748, 19265, 27028,
                                                                       13625, 13850, 20735,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 28288, 0, 3,
                                                                       27028, 19475, 27308,
                                                                       13850, 14075, 21050,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 28708, 0, 3,
                                                                       27308, 19685, 27588,
                                                                       14075, 14300, 21365,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 29128, 0, 3,
                                                                       27868, 20735, 28288,
                                                                       14750, 15065, 22562,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 29716, 0, 3,
                                                                       28288, 21050, 28708,
                                                                       15065, 15380, 23003,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 30304, 0, 3,
                                                                       29128, 22562, 29716,
                                                                       16010, 16430, 24620,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 31088, 3, 17270,
                                                                       17291, 25208, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 31124, 3, 17291,
                                                                       17312, 25236, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 31160, 3, 17312,
                                                                       17333, 25264, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 31196, 3, 17333,
                                                                       17354, 25292, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 31232, 3, 17354,
                                                                       17375, 25320, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 31268, 3, 17375,
                                                                       17396, 25348, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 31304, 3, 17396,
                                                                       17417, 25376, ncols,
                                                                       gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 31340, 0, 3,
                                                                       31088, 25208, 31124,
                                                                       17459, 17522, 25404,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 31448, 0, 3,
                                                                       31124, 25236, 31160,
                                                                       17522, 17585, 25488,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 31556, 0, 3,
                                                                       31160, 25264, 31196,
                                                                       17585, 17648, 25572,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 31664, 0, 3,
                                                                       31196, 25292, 31232,
                                                                       17648, 17711, 25656,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 31772, 0, 3,
                                                                       31232, 25320, 31268,
                                                                       17711, 17774, 25740,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 31880, 0, 3,
                                                                       31268, 25348, 31304,
                                                                       17774, 17837, 25824,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 31988, 0, 3,
                                                                       31340, 25404, 31448,
                                                                       17963, 18089, 25908,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 32204, 0, 3,
                                                                       31448, 25488, 31556,
                                                                       18089, 18215, 26076,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 32420, 0, 3,
                                                                       31556, 25572, 31664,
                                                                       18215, 18341, 26244,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 32636, 0, 3,
                                                                       31664, 25656, 31772,
                                                                       18341, 18467, 26412,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 32852, 0, 3,
                                                                       31772, 25740, 31880,
                                                                       18467, 18593, 26580,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 33068, 0, 3,
                                                                       31988, 25908, 32204,
                                                                       18845, 19055, 26748,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 33428, 0, 3,
                                                                       32204, 26076, 32420,
                                                                       19055, 19265, 27028,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 33788, 0, 3,
                                                                       32420, 26244, 32636,
                                                                       19265, 19475, 27308,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 34148, 0, 3,
                                                                       32636, 26412, 32852,
                                                                       19475, 19685, 27588,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 34508, 0, 3,
                                                                       33068, 26748, 33428,
                                                                       20105, 20420, 27868,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 35048, 0, 3,
                                                                       33428, 27028, 33788,
                                                                       20420, 20735, 28288,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 35588, 0, 3,
                                                                       33788, 27308, 34148,
                                                                       20735, 21050, 28708,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 36128, 0, 3,
                                                                       34508, 27868, 35048,
                                                                       21680, 22121, 29128,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 36884, 0, 3,
                                                                       35048, 28288, 35588,
                                                                       22121, 22562, 29716,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 37640, 0, 3,
                                                                       36128, 29128, 36884,
                                                                       23444, 24032, 30304,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 38648, 36128, 756, ncols);

                    simdfunc::contract_primitives(buffer, 39719, 37640, 1008, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 39404, 38648, 21, 1, nmax);

        simdtrf::transform_k_inner(buffer, 40727, 39719, 28, 1, nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 41147, 39404, 40727, 15,
                                             nmax);

        simdtrf::transform_p_inner(buffer, 42092, 41147, 21, 15, nmax);

        simdtrf::transform_h_outer(values + n * npairs, nvalues, buffer, 42092, 45, nmax);
    }

    for (size_t m = 0; m < 495; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
