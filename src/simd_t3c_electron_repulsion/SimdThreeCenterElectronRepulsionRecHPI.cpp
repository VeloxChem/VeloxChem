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


#include "SimdThreeCenterElectronRepulsionRecHPI.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferHP.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformI.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_hpi_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_hpi_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 29023, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 429 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 29023, 25376, 1645, dimensions);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 722, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 725, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 728, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 731, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 734, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 737, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 740, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 743, 3, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 746, 3, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 749, 3, 19, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 752, 3, 20, ncols,
                                                                       p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 755, 3, 10, 27,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 764, 3, 11, 30,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 773, 3, 12, 33,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 782, 3, 13, 36,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 791, 3, 14, 39,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 800, 3, 15, 42,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 809, 3, 16, 45,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 818, 3, 17, 48,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 827, 3, 18, 51,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 836, 3, 19, 54,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 845, 3, 27, 69,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 863, 3, 30, 75,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 881, 3, 33, 81,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 899, 3, 36, 87,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 917, 3, 39, 93,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 935, 3, 42, 99,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 953, 3, 45, 105,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 971, 3, 48, 111,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 989, 3, 51, 117,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1007, 3, 69, 143,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1037, 3, 75, 153,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1067, 3, 81, 163,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1097, 3, 87, 173,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1127, 3, 93, 183,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1157, 3, 99, 193,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1187, 3, 105, 203,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1217, 3, 111, 213,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1247, 3, 143, 253,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1292, 3, 153, 268,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1337, 3, 163, 283,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1382, 3, 173, 298,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1427, 3, 183, 313,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1472, 3, 193, 328,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1517, 3, 203, 343,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1562, 3, 253, 400,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1625, 3, 268, 421,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1688, 3, 283, 442,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1751, 3, 298, 463,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1814, 3, 313, 484,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1877, 3, 328, 505,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 1940, 3, 400, 582,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2024, 3, 421, 610,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2108, 3, 442, 638,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2192, 3, 463, 666,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2276, 3, 484, 694,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2360, 3, 8, 9,
                                                                       722, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2366, 3, 9, 10,
                                                                       725, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2372, 3, 10, 11,
                                                                       728, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2378, 3, 11, 12,
                                                                       731, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2384, 3, 12, 13,
                                                                       734, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2390, 3, 13, 14,
                                                                       737, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2396, 3, 14, 15,
                                                                       740, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2402, 3, 15, 16,
                                                                       743, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2408, 3, 16, 17,
                                                                       746, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2414, 3, 17, 18,
                                                                       749, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2420, 3, 18, 19,
                                                                       752, ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2426, 0, 3, 2360,
                                                                       722, 2366, 755, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2444, 0, 3, 2366,
                                                                       725, 2372, 764, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2462, 0, 3, 2372,
                                                                       728, 2378, 773, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2480, 0, 3, 2378,
                                                                       731, 2384, 782, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2498, 0, 3, 2384,
                                                                       734, 2390, 791, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2516, 0, 3, 2390,
                                                                       737, 2396, 800, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2534, 0, 3, 2396,
                                                                       740, 2402, 809, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2552, 0, 3, 2402,
                                                                       743, 2408, 818, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2570, 0, 3, 2408,
                                                                       746, 2414, 827, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2588, 0, 3, 2414,
                                                                       749, 2420, 836, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2606, 0, 3, 2426,
                                                                       755, 2444, 57, 63, 845,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2642, 0, 3, 2444,
                                                                       764, 2462, 63, 69, 863,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2678, 0, 3, 2462,
                                                                       773, 2480, 69, 75, 881,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2714, 0, 3, 2480,
                                                                       782, 2498, 75, 81, 899,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2750, 0, 3, 2498,
                                                                       791, 2516, 81, 87, 917,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2786, 0, 3, 2516,
                                                                       800, 2534, 87, 93, 935,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2822, 0, 3, 2534,
                                                                       809, 2552, 93, 99, 953,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2858, 0, 3, 2552,
                                                                       818, 2570, 99, 105, 971,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2894, 0, 3, 2570,
                                                                       827, 2588, 105, 111, 989,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2930, 0, 3, 2606,
                                                                       845, 2642, 123, 133, 1007,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2990, 0, 3, 2642,
                                                                       863, 2678, 133, 143, 1037,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3050, 0, 3, 2678,
                                                                       881, 2714, 143, 153, 1067,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3110, 0, 3, 2714,
                                                                       899, 2750, 153, 163, 1097,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3170, 0, 3, 2750,
                                                                       917, 2786, 163, 173, 1127,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3230, 0, 3, 2786,
                                                                       935, 2822, 173, 183, 1157,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3290, 0, 3, 2822,
                                                                       953, 2858, 183, 193, 1187,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3350, 0, 3, 2858,
                                                                       971, 2894, 193, 203, 1217,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 3410, 0, 3, 2930,
                                                                       1007, 2990, 223, 238,
                                                                       1247, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 3500, 0, 3, 2990,
                                                                       1037, 3050, 238, 253,
                                                                       1292, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 3590, 0, 3, 3050,
                                                                       1067, 3110, 253, 268,
                                                                       1337, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 3680, 0, 3, 3110,
                                                                       1097, 3170, 268, 283,
                                                                       1382, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 3770, 0, 3, 3170,
                                                                       1127, 3230, 283, 298,
                                                                       1427, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 3860, 0, 3, 3230,
                                                                       1157, 3290, 298, 313,
                                                                       1472, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 3950, 0, 3, 3290,
                                                                       1187, 3350, 313, 328,
                                                                       1517, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 4040, 0, 3, 3410,
                                                                       1247, 3500, 358, 379,
                                                                       1562, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 4166, 0, 3, 3500,
                                                                       1292, 3590, 379, 400,
                                                                       1625, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 4292, 0, 3, 3590,
                                                                       1337, 3680, 400, 421,
                                                                       1688, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 4418, 0, 3, 3680,
                                                                       1382, 3770, 421, 442,
                                                                       1751, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 4544, 0, 3, 3770,
                                                                       1427, 3860, 442, 463,
                                                                       1814, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 4670, 0, 3, 3860,
                                                                       1472, 3950, 463, 484,
                                                                       1877, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 4796, 0, 3, 4040,
                                                                       1562, 4166, 526, 554,
                                                                       1940, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 4964, 0, 3, 4166,
                                                                       1625, 4292, 554, 582,
                                                                       2024, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 5132, 0, 3, 4292,
                                                                       1688, 4418, 582, 610,
                                                                       2108, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 5300, 0, 3, 4418,
                                                                       1751, 4544, 610, 638,
                                                                       2192, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 5468, 0, 3, 4544,
                                                                       1814, 4670, 638, 666,
                                                                       2276, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5636, 3, 722, 725,
                                                                       2372, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5646, 3, 725, 728,
                                                                       2378, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5656, 3, 728, 731,
                                                                       2384, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5666, 3, 731, 734,
                                                                       2390, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5676, 3, 734, 737,
                                                                       2396, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5686, 3, 737, 740,
                                                                       2402, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5696, 3, 740, 743,
                                                                       2408, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5706, 3, 743, 746,
                                                                       2414, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5716, 3, 746, 749,
                                                                       2420, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5726, 0, 3, 5636,
                                                                       2372, 5646, 755, 764,
                                                                       2462, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5756, 0, 3, 5646,
                                                                       2378, 5656, 764, 773,
                                                                       2480, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5786, 0, 3, 5656,
                                                                       2384, 5666, 773, 782,
                                                                       2498, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5816, 0, 3, 5666,
                                                                       2390, 5676, 782, 791,
                                                                       2516, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5846, 0, 3, 5676,
                                                                       2396, 5686, 791, 800,
                                                                       2534, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5876, 0, 3, 5686,
                                                                       2402, 5696, 800, 809,
                                                                       2552, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5906, 0, 3, 5696,
                                                                       2408, 5706, 809, 818,
                                                                       2570, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5936, 0, 3, 5706,
                                                                       2414, 5716, 818, 827,
                                                                       2588, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 5966, 0, 3, 5726,
                                                                       2462, 5756, 845, 863,
                                                                       2678, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6026, 0, 3, 5756,
                                                                       2480, 5786, 863, 881,
                                                                       2714, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6086, 0, 3, 5786,
                                                                       2498, 5816, 881, 899,
                                                                       2750, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6146, 0, 3, 5816,
                                                                       2516, 5846, 899, 917,
                                                                       2786, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6206, 0, 3, 5846,
                                                                       2534, 5876, 917, 935,
                                                                       2822, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6266, 0, 3, 5876,
                                                                       2552, 5906, 935, 953,
                                                                       2858, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6326, 0, 3, 5906,
                                                                       2570, 5936, 953, 971,
                                                                       2894, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 6386, 0, 3, 5966,
                                                                       2678, 6026, 1007, 1037,
                                                                       3050, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 6486, 0, 3, 6026,
                                                                       2714, 6086, 1037, 1067,
                                                                       3110, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 6586, 0, 3, 6086,
                                                                       2750, 6146, 1067, 1097,
                                                                       3170, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 6686, 0, 3, 6146,
                                                                       2786, 6206, 1097, 1127,
                                                                       3230, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 6786, 0, 3, 6206,
                                                                       2822, 6266, 1127, 1157,
                                                                       3290, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 6886, 0, 3, 6266,
                                                                       2858, 6326, 1157, 1187,
                                                                       3350, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 6986, 0, 3, 6386,
                                                                       3050, 6486, 1247, 1292,
                                                                       3590, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 7136, 0, 3, 6486,
                                                                       3110, 6586, 1292, 1337,
                                                                       3680, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 7286, 0, 3, 6586,
                                                                       3170, 6686, 1337, 1382,
                                                                       3770, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 7436, 0, 3, 6686,
                                                                       3230, 6786, 1382, 1427,
                                                                       3860, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 7586, 0, 3, 6786,
                                                                       3290, 6886, 1427, 1472,
                                                                       3950, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 7736, 0, 3, 6986,
                                                                       3590, 7136, 1562, 1625,
                                                                       4292, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 7946, 0, 3, 7136,
                                                                       3680, 7286, 1625, 1688,
                                                                       4418, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 8156, 0, 3, 7286,
                                                                       3770, 7436, 1688, 1751,
                                                                       4544, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 8366, 0, 3, 7436,
                                                                       3860, 7586, 1751, 1814,
                                                                       4670, ncols, gamma, p,
                                                                       q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 8576, 0, 3, 7736,
                                                                       4292, 7946, 1940, 2024,
                                                                       5132, ncols, gamma, p,
                                                                       q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 8856, 0, 3, 7946,
                                                                       4418, 8156, 2024, 2108,
                                                                       5300, ncols, gamma, p,
                                                                       q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 9136, 0, 3, 8156,
                                                                       4544, 8366, 2108, 2192,
                                                                       5468, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9416, 3, 2360,
                                                                       2366, 5636, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9431, 3, 2366,
                                                                       2372, 5646, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9446, 3, 2372,
                                                                       2378, 5656, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9461, 3, 2378,
                                                                       2384, 5666, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9476, 3, 2384,
                                                                       2390, 5676, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9491, 3, 2390,
                                                                       2396, 5686, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9506, 3, 2396,
                                                                       2402, 5696, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9521, 3, 2402,
                                                                       2408, 5706, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9536, 3, 2408,
                                                                       2414, 5716, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 9551, 0, 3, 9416,
                                                                       5636, 9431, 2426, 2444,
                                                                       5726, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 9596, 0, 3, 9431,
                                                                       5646, 9446, 2444, 2462,
                                                                       5756, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 9641, 0, 3, 9446,
                                                                       5656, 9461, 2462, 2480,
                                                                       5786, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 9686, 0, 3, 9461,
                                                                       5666, 9476, 2480, 2498,
                                                                       5816, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 9731, 0, 3, 9476,
                                                                       5676, 9491, 2498, 2516,
                                                                       5846, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 9776, 0, 3, 9491,
                                                                       5686, 9506, 2516, 2534,
                                                                       5876, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 9821, 0, 3, 9506,
                                                                       5696, 9521, 2534, 2552,
                                                                       5906, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 9866, 0, 3, 9521,
                                                                       5706, 9536, 2552, 2570,
                                                                       5936, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 9911, 0, 3, 9551,
                                                                       5726, 9596, 2606, 2642,
                                                                       5966, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 10001, 0, 3, 9596,
                                                                       5756, 9641, 2642, 2678,
                                                                       6026, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 10091, 0, 3, 9641,
                                                                       5786, 9686, 2678, 2714,
                                                                       6086, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 10181, 0, 3, 9686,
                                                                       5816, 9731, 2714, 2750,
                                                                       6146, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 10271, 0, 3, 9731,
                                                                       5846, 9776, 2750, 2786,
                                                                       6206, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 10361, 0, 3, 9776,
                                                                       5876, 9821, 2786, 2822,
                                                                       6266, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 10451, 0, 3, 9821,
                                                                       5906, 9866, 2822, 2858,
                                                                       6326, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 10541, 0, 3, 9911,
                                                                       5966, 10001, 2930, 2990,
                                                                       6386, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 10691, 0, 3,
                                                                       10001, 6026, 10091, 2990,
                                                                       3050, 6486, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 10841, 0, 3,
                                                                       10091, 6086, 10181, 3050,
                                                                       3110, 6586, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 10991, 0, 3,
                                                                       10181, 6146, 10271, 3110,
                                                                       3170, 6686, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 11141, 0, 3,
                                                                       10271, 6206, 10361, 3170,
                                                                       3230, 6786, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 11291, 0, 3,
                                                                       10361, 6266, 10451, 3230,
                                                                       3290, 6886, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 11441, 0, 3,
                                                                       10541, 6386, 10691, 3410,
                                                                       3500, 6986, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 11666, 0, 3,
                                                                       10691, 6486, 10841, 3500,
                                                                       3590, 7136, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 11891, 0, 3,
                                                                       10841, 6586, 10991, 3590,
                                                                       3680, 7286, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 12116, 0, 3,
                                                                       10991, 6686, 11141, 3680,
                                                                       3770, 7436, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 12341, 0, 3,
                                                                       11141, 6786, 11291, 3770,
                                                                       3860, 7586, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 12566, 0, 3,
                                                                       11441, 6986, 11666, 4040,
                                                                       4166, 7736, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 12881, 0, 3,
                                                                       11666, 7136, 11891, 4166,
                                                                       4292, 7946, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 13196, 0, 3,
                                                                       11891, 7286, 12116, 4292,
                                                                       4418, 8156, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 13511, 0, 3,
                                                                       12116, 7436, 12341, 4418,
                                                                       4544, 8366, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 13826, 0, 3,
                                                                       12566, 7736, 12881, 4796,
                                                                       4964, 8576, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 14246, 0, 3,
                                                                       12881, 7946, 13196, 4964,
                                                                       5132, 8856, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 14666, 0, 3,
                                                                       13196, 8156, 13511, 5132,
                                                                       5300, 9136, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 15086, 3, 5636,
                                                                       5646, 9446, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 15107, 3, 5646,
                                                                       5656, 9461, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 15128, 3, 5656,
                                                                       5666, 9476, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 15149, 3, 5666,
                                                                       5676, 9491, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 15170, 3, 5676,
                                                                       5686, 9506, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 15191, 3, 5686,
                                                                       5696, 9521, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 15212, 3, 5696,
                                                                       5706, 9536, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 15233, 0, 3,
                                                                       15086, 9446, 15107, 5726,
                                                                       5756, 9641, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 15296, 0, 3,
                                                                       15107, 9461, 15128, 5756,
                                                                       5786, 9686, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 15359, 0, 3,
                                                                       15128, 9476, 15149, 5786,
                                                                       5816, 9731, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 15422, 0, 3,
                                                                       15149, 9491, 15170, 5816,
                                                                       5846, 9776, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 15485, 0, 3,
                                                                       15170, 9506, 15191, 5846,
                                                                       5876, 9821, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 15548, 0, 3,
                                                                       15191, 9521, 15212, 5876,
                                                                       5906, 9866, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 15611, 0, 3,
                                                                       15233, 9641, 15296, 5966,
                                                                       6026, 10091, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 15737, 0, 3,
                                                                       15296, 9686, 15359, 6026,
                                                                       6086, 10181, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 15863, 0, 3,
                                                                       15359, 9731, 15422, 6086,
                                                                       6146, 10271, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 15989, 0, 3,
                                                                       15422, 9776, 15485, 6146,
                                                                       6206, 10361, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 16115, 0, 3,
                                                                       15485, 9821, 15548, 6206,
                                                                       6266, 10451, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 16241, 0, 3,
                                                                       15611, 10091, 15737, 6386,
                                                                       6486, 10841, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 16451, 0, 3,
                                                                       15737, 10181, 15863, 6486,
                                                                       6586, 10991, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 16661, 0, 3,
                                                                       15863, 10271, 15989, 6586,
                                                                       6686, 11141, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 16871, 0, 3,
                                                                       15989, 10361, 16115, 6686,
                                                                       6786, 11291, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 17081, 0, 3,
                                                                       16241, 10841, 16451, 6986,
                                                                       7136, 11891, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 17396, 0, 3,
                                                                       16451, 10991, 16661, 7136,
                                                                       7286, 12116, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 17711, 0, 3,
                                                                       16661, 11141, 16871, 7286,
                                                                       7436, 12341, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 18026, 0, 3,
                                                                       17081, 11891, 17396, 7736,
                                                                       7946, 13196, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 18467, 0, 3,
                                                                       17396, 12116, 17711, 7946,
                                                                       8156, 13511, ncols, gamma,
                                                                       p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 18908, 0, 3,
                                                                       18026, 13196, 18467, 8576,
                                                                       8856, 14666, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 19496, 3, 9416,
                                                                       9431, 15086, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 19524, 3, 9431,
                                                                       9446, 15107, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 19552, 3, 9446,
                                                                       9461, 15128, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 19580, 3, 9461,
                                                                       9476, 15149, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 19608, 3, 9476,
                                                                       9491, 15170, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 19636, 3, 9491,
                                                                       9506, 15191, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 19664, 3, 9506,
                                                                       9521, 15212, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 19692, 0, 3,
                                                                       19496, 15086, 19524, 9551,
                                                                       9596, 15233, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 19776, 0, 3,
                                                                       19524, 15107, 19552, 9596,
                                                                       9641, 15296, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 19860, 0, 3,
                                                                       19552, 15128, 19580, 9641,
                                                                       9686, 15359, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 19944, 0, 3,
                                                                       19580, 15149, 19608, 9686,
                                                                       9731, 15422, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 20028, 0, 3,
                                                                       19608, 15170, 19636, 9731,
                                                                       9776, 15485, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 20112, 0, 3,
                                                                       19636, 15191, 19664, 9776,
                                                                       9821, 15548, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 20196, 0, 3,
                                                                       19692, 15233, 19776, 9911,
                                                                       10001, 15611, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 20364, 0, 3,
                                                                       19776, 15296, 19860,
                                                                       10001, 10091, 15737,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 20532, 0, 3,
                                                                       19860, 15359, 19944,
                                                                       10091, 10181, 15863,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 20700, 0, 3,
                                                                       19944, 15422, 20028,
                                                                       10181, 10271, 15989,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 20868, 0, 3,
                                                                       20028, 15485, 20112,
                                                                       10271, 10361, 16115,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 21036, 0, 3,
                                                                       20196, 15611, 20364,
                                                                       10541, 10691, 16241,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 21316, 0, 3,
                                                                       20364, 15737, 20532,
                                                                       10691, 10841, 16451,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 21596, 0, 3,
                                                                       20532, 15863, 20700,
                                                                       10841, 10991, 16661,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 21876, 0, 3,
                                                                       20700, 15989, 20868,
                                                                       10991, 11141, 16871,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 22156, 0, 3,
                                                                       21036, 16241, 21316,
                                                                       11441, 11666, 17081,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 22576, 0, 3,
                                                                       21316, 16451, 21596,
                                                                       11666, 11891, 17396,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 22996, 0, 3,
                                                                       21596, 16661, 21876,
                                                                       11891, 12116, 17711,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 23416, 0, 3,
                                                                       22156, 17081, 22576,
                                                                       12566, 12881, 18026,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 24004, 0, 3,
                                                                       22576, 17396, 22996,
                                                                       12881, 13196, 18467,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 24592, 0, 3,
                                                                       23416, 18026, 24004,
                                                                       13826, 14246, 18908,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 25376, 23416, 588, ncols);

                    simdfunc::contract_primitives(buffer, 26237, 24592, 784, ncols);
                }
            }
        }

        simdtrf::transform_i_inner(buffer, 25964, 25376, 21, 1, nmax);

        simdtrf::transform_i_inner(buffer, 27021, 26237, 28, 1, nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 27385, 25964, 27021, 13,
                                             nmax);

        simdtrf::transform_p_inner(buffer, 28204, 27385, 21, 13, nmax);

        simdtrf::transform_h_outer(values + n * npairs, nvalues, buffer, 28204, 39, nmax);
    }

    for (size_t m = 0; m < 429; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
