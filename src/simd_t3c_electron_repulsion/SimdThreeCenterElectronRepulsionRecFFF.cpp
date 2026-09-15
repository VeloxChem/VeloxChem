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


#include "SimdThreeCenterElectronRepulsionRecFFF.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferFD.hpp"
#include "SimdTransferFF.hpp"
#include "SimdTransferFP.hpp"
#include "SimdTransferGD.hpp"
#include "SimdTransferGP.hpp"
#include "SimdTransferHP.hpp"
#include "SimdTransformF.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_fff_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_fff_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 9344, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 343 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 9344, 4880, 1062, dimensions);

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
                                                        5, 6, 7, 8, 9}, ncols, fj, 6, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 17, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 20, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 23, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 26, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 29, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 32, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 35, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 38, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 41, 0, 3, 8, 9,
                                                                       17, 20, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 47, 0, 3, 9, 10,
                                                                       20, 23, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 53, 0, 3, 10, 11,
                                                                       23, 26, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 59, 0, 3, 11, 12,
                                                                       26, 29, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 65, 0, 3, 12, 13,
                                                                       29, 32, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 71, 0, 3, 13, 14,
                                                                       32, 35, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 77, 0, 3, 14, 15,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 83, 0, 3, 17, 20,
                                                                       41, 47, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 93, 0, 3, 20, 23,
                                                                       47, 53, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 103, 0, 3, 23, 26,
                                                                       53, 59, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 113, 0, 3, 26, 29,
                                                                       59, 65, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 123, 0, 3, 29, 32,
                                                                       65, 71, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 133, 0, 3, 32, 35,
                                                                       71, 77, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 143, 0, 3, 41, 47,
                                                                       83, 93, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 158, 0, 3, 47, 53,
                                                                       93, 103, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 173, 0, 3, 53, 59,
                                                                       103, 113, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 188, 0, 3, 59, 65,
                                                                       113, 123, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 203, 0, 3, 65, 71,
                                                                       123, 133, ncols, gamma, p,
                                                                       q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 218, 0, 3, 83, 93,
                                                                       143, 158, ncols, gamma, p,
                                                                       q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 239, 0, 3, 93,
                                                                       103, 158, 173, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 260, 0, 3, 103,
                                                                       113, 173, 188, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 281, 0, 3, 113,
                                                                       123, 188, 203, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 302, 0, 3, 143,
                                                                       158, 218, 239, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 330, 0, 3, 158,
                                                                       173, 239, 260, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 358, 0, 3, 173,
                                                                       188, 260, 281, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 386, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 389, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 392, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 395, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 398, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 401, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 404, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 407, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 410, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 413, 3, 8, 17,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 422, 3, 9, 20,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 431, 3, 10, 23,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 440, 3, 11, 26,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 449, 3, 12, 29,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 458, 3, 13, 32,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 467, 3, 14, 35,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 476, 3, 15, 38,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 485, 3, 17, 41,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 503, 3, 20, 47,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 521, 3, 23, 53,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 539, 3, 26, 59,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 557, 3, 29, 65,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 575, 3, 32, 71,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 593, 3, 35, 77,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 611, 3, 41, 83,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 641, 3, 47, 93,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 671, 3, 53, 103,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 701, 3, 59, 113,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 731, 3, 65, 123,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 761, 3, 71, 133,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 791, 3, 83, 143,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 836, 3, 93, 158,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 881, 3, 103, 173,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 926, 3, 113, 188,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 971, 3, 123, 203,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1016, 3, 143, 218,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1079, 3, 158, 239,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1142, 3, 173, 260,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1205, 3, 188, 281,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 1268, 3, 218, 302,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 1352, 3, 239, 330,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 1436, 3, 260, 358,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1520, 3, 8, 9,
                                                                       392, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1526, 3, 9, 10,
                                                                       395, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1532, 3, 10, 11,
                                                                       398, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1538, 3, 11, 12,
                                                                       401, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1544, 3, 12, 13,
                                                                       404, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1550, 3, 13, 14,
                                                                       407, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1556, 3, 14, 15,
                                                                       410, ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1562, 0, 3, 1520,
                                                                       392, 1526, 431, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1580, 0, 3, 1526,
                                                                       395, 1532, 440, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1598, 0, 3, 1532,
                                                                       398, 1538, 449, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1616, 0, 3, 1538,
                                                                       401, 1544, 458, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1634, 0, 3, 1544,
                                                                       404, 1550, 467, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1652, 0, 3, 1550,
                                                                       407, 1556, 476, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1670, 0, 3, 1562,
                                                                       431, 1580, 41, 47, 521,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1706, 0, 3, 1580,
                                                                       440, 1598, 47, 53, 539,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1742, 0, 3, 1598,
                                                                       449, 1616, 53, 59, 557,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1778, 0, 3, 1616,
                                                                       458, 1634, 59, 65, 575,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1814, 0, 3, 1634,
                                                                       467, 1652, 65, 71, 593,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 1850, 0, 3, 1670,
                                                                       521, 1706, 83, 93, 671,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 1910, 0, 3, 1706,
                                                                       539, 1742, 93, 103, 701,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 1970, 0, 3, 1742,
                                                                       557, 1778, 103, 113, 731,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2030, 0, 3, 1778,
                                                                       575, 1814, 113, 123, 761,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 2090, 0, 3, 1850,
                                                                       671, 1910, 143, 158, 881,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 2180, 0, 3, 1910,
                                                                       701, 1970, 158, 173, 926,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 2270, 0, 3, 1970,
                                                                       731, 2030, 173, 188, 971,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 2360, 0, 3, 2090,
                                                                       881, 2180, 218, 239, 1142,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 2486, 0, 3, 2180,
                                                                       926, 2270, 239, 260, 1205,
                                                                       ncols, gamma, p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 2612, 0, 3, 2360,
                                                                       1142, 2486, 302, 330,
                                                                       1436, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2780, 3, 386, 389,
                                                                       1520, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2790, 3, 389, 392,
                                                                       1526, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2800, 3, 392, 395,
                                                                       1532, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2810, 3, 395, 398,
                                                                       1538, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2820, 3, 398, 401,
                                                                       1544, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2830, 3, 401, 404,
                                                                       1550, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2840, 3, 404, 407,
                                                                       1556, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2850, 0, 3, 2780,
                                                                       1520, 2790, 413, 422,
                                                                       1562, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2880, 0, 3, 2790,
                                                                       1526, 2800, 422, 431,
                                                                       1580, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2910, 0, 3, 2800,
                                                                       1532, 2810, 431, 440,
                                                                       1598, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2940, 0, 3, 2810,
                                                                       1538, 2820, 440, 449,
                                                                       1616, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2970, 0, 3, 2820,
                                                                       1544, 2830, 449, 458,
                                                                       1634, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 3000, 0, 3, 2830,
                                                                       1550, 2840, 458, 467,
                                                                       1652, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 3030, 0, 3, 2850,
                                                                       1562, 2880, 485, 503,
                                                                       1670, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 3090, 0, 3, 2880,
                                                                       1580, 2910, 503, 521,
                                                                       1706, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 3150, 0, 3, 2910,
                                                                       1598, 2940, 521, 539,
                                                                       1742, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 3210, 0, 3, 2940,
                                                                       1616, 2970, 539, 557,
                                                                       1778, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 3270, 0, 3, 2970,
                                                                       1634, 3000, 557, 575,
                                                                       1814, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 3330, 0, 3, 3030,
                                                                       1670, 3090, 611, 641,
                                                                       1850, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 3430, 0, 3, 3090,
                                                                       1706, 3150, 641, 671,
                                                                       1910, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 3530, 0, 3, 3150,
                                                                       1742, 3210, 671, 701,
                                                                       1970, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 3630, 0, 3, 3210,
                                                                       1778, 3270, 701, 731,
                                                                       2030, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 3730, 0, 3, 3330,
                                                                       1850, 3430, 791, 836,
                                                                       2090, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 3880, 0, 3, 3430,
                                                                       1910, 3530, 836, 881,
                                                                       2180, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 4030, 0, 3, 3530,
                                                                       1970, 3630, 881, 926,
                                                                       2270, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 4180, 0, 3, 3730,
                                                                       2090, 3880, 1016, 1079,
                                                                       2360, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 4390, 0, 3, 3880,
                                                                       2180, 4030, 1079, 1142,
                                                                       2486, ncols, gamma, p,
                                                                       q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 4600, 0, 3, 4180,
                                                                       2360, 4390, 1268, 1352,
                                                                       2612, ncols, gamma, p,
                                                                       q);

                    simdfunc::contract_primitives(buffer, 4880, 3330, 100, ncols);

                    simdfunc::contract_primitives(buffer, 5050, 3730, 150, ncols);

                    simdfunc::contract_primitives(buffer, 5305, 4180, 210, ncols);

                    simdfunc::contract_primitives(buffer, 5662, 4600, 280, ncols);
                }
            }
        }

        simdtrf::transform_f_inner(buffer, 4980, 4880, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 5200, 5050, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 5515, 5305, 21, 1, nmax);

        simdtrf::transform_f_inner(buffer, 5942, 5662, 28, 1, nmax);

        simdtrf::compute_hrr_fp_out_of_first(buffer, coordinates, 6138, 4980, 5200, 7, nmax);

        simdtrf::compute_hrr_gp_out_of_first(buffer, coordinates, 6348, 5200, 5515, 7, nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 6663, 5515, 5942, 7, nmax);

        simdtrf::compute_hrr_fd_out_of_first(buffer, coordinates, 7104, 6138, 6348, 7, nmax);

        simdtrf::compute_hrr_gd_out_of_first(buffer, coordinates, 7524, 6348, 6663, 7, nmax);

        simdtrf::compute_hrr_ff_out_of_first(buffer, coordinates, 8154, 7104, 7524, 7, nmax);

        simdtrf::transform_f_inner(buffer, 8854, 8154, 10, 7, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 8854, 49, nmax);
    }

    for (size_t m = 0; m < 343; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
