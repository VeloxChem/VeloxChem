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


#include "SimdThreeCenterElectronRepulsionRecGPI.hpp"

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
#include "SimdTransferGP.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformI.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_gpi_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_gpi_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 18390, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 351 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 18390, 15744, 1203, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 7, 3, 11,
                                                             ncols, fj, 6, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 20, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 23, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 26, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 29, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 32, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 35, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 38, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 41, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 44, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 47, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 50, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 53, 0, 3, 8, 9,
                                                                       20, 23, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 59, 0, 3, 9, 10,
                                                                       23, 26, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 65, 0, 3, 10, 11,
                                                                       26, 29, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 71, 0, 3, 11, 12,
                                                                       29, 32, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 77, 0, 3, 12, 13,
                                                                       32, 35, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 83, 0, 3, 13, 14,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 89, 0, 3, 14, 15,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 95, 0, 3, 15, 16,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 101, 0, 3, 16, 17,
                                                                       44, 47, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 107, 0, 3, 17, 18,
                                                                       47, 50, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 113, 0, 3, 20, 23,
                                                                       53, 59, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 123, 0, 3, 23, 26,
                                                                       59, 65, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 133, 0, 3, 26, 29,
                                                                       65, 71, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 143, 0, 3, 29, 32,
                                                                       71, 77, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 153, 0, 3, 32, 35,
                                                                       77, 83, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 163, 0, 3, 35, 38,
                                                                       83, 89, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 173, 0, 3, 38, 41,
                                                                       89, 95, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 183, 0, 3, 41, 44,
                                                                       95, 101, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 193, 0, 3, 44, 47,
                                                                       101, 107, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 203, 0, 3, 53, 59,
                                                                       113, 123, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 218, 0, 3, 59, 65,
                                                                       123, 133, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 233, 0, 3, 65, 71,
                                                                       133, 143, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 248, 0, 3, 71, 77,
                                                                       143, 153, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 263, 0, 3, 77, 83,
                                                                       153, 163, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 278, 0, 3, 83, 89,
                                                                       163, 173, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 293, 0, 3, 89, 95,
                                                                       173, 183, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 308, 0, 3, 95,
                                                                       101, 183, 193, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 323, 0, 3, 113,
                                                                       123, 203, 218, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 344, 0, 3, 123,
                                                                       133, 218, 233, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 365, 0, 3, 133,
                                                                       143, 233, 248, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 386, 0, 3, 143,
                                                                       153, 248, 263, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 407, 0, 3, 153,
                                                                       163, 263, 278, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 428, 0, 3, 163,
                                                                       173, 278, 293, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 449, 0, 3, 173,
                                                                       183, 293, 308, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 470, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 473, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 476, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 479, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 482, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 485, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 488, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 491, 3, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 494, 3, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 497, 3, 19, ncols,
                                                                       p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 500, 3, 10, 26,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 509, 3, 11, 29,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 518, 3, 12, 32,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 527, 3, 13, 35,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 536, 3, 14, 38,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 545, 3, 15, 41,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 554, 3, 16, 44,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 563, 3, 17, 47,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 572, 3, 18, 50,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 581, 3, 26, 65,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 599, 3, 29, 71,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 617, 3, 32, 77,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 635, 3, 35, 83,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 653, 3, 38, 89,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 671, 3, 41, 95,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 689, 3, 44, 101,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 707, 3, 47, 107,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 725, 3, 65, 133,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 755, 3, 71, 143,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 785, 3, 77, 153,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 815, 3, 83, 163,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 845, 3, 89, 173,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 875, 3, 95, 183,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 905, 3, 101, 193,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 935, 3, 133, 233,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 980, 3, 143, 248,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1025, 3, 153, 263,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1070, 3, 163, 278,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1115, 3, 173, 293,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1160, 3, 183, 308,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1205, 3, 233, 365,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1268, 3, 248, 386,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1331, 3, 263, 407,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1394, 3, 278, 428,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1457, 3, 293, 449,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1520, 3, 8, 9,
                                                                       470, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1526, 3, 9, 10,
                                                                       473, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1532, 3, 10, 11,
                                                                       476, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1538, 3, 11, 12,
                                                                       479, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1544, 3, 12, 13,
                                                                       482, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1550, 3, 13, 14,
                                                                       485, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1556, 3, 14, 15,
                                                                       488, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1562, 3, 15, 16,
                                                                       491, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1568, 3, 16, 17,
                                                                       494, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1574, 3, 17, 18,
                                                                       497, ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1580, 0, 3, 1520,
                                                                       470, 1526, 500, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1598, 0, 3, 1526,
                                                                       473, 1532, 509, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1616, 0, 3, 1532,
                                                                       476, 1538, 518, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1634, 0, 3, 1538,
                                                                       479, 1544, 527, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1652, 0, 3, 1544,
                                                                       482, 1550, 536, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1670, 0, 3, 1550,
                                                                       485, 1556, 545, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1688, 0, 3, 1556,
                                                                       488, 1562, 554, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1706, 0, 3, 1562,
                                                                       491, 1568, 563, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1724, 0, 3, 1568,
                                                                       494, 1574, 572, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1742, 0, 3, 1580,
                                                                       500, 1598, 53, 59, 581,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1778, 0, 3, 1598,
                                                                       509, 1616, 59, 65, 599,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1814, 0, 3, 1616,
                                                                       518, 1634, 65, 71, 617,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1850, 0, 3, 1634,
                                                                       527, 1652, 71, 77, 635,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1886, 0, 3, 1652,
                                                                       536, 1670, 77, 83, 653,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1922, 0, 3, 1670,
                                                                       545, 1688, 83, 89, 671,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1958, 0, 3, 1688,
                                                                       554, 1706, 89, 95, 689,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1994, 0, 3, 1706,
                                                                       563, 1724, 95, 101, 707,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2030, 0, 3, 1742,
                                                                       581, 1778, 113, 123, 725,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2090, 0, 3, 1778,
                                                                       599, 1814, 123, 133, 755,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2150, 0, 3, 1814,
                                                                       617, 1850, 133, 143, 785,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2210, 0, 3, 1850,
                                                                       635, 1886, 143, 153, 815,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2270, 0, 3, 1886,
                                                                       653, 1922, 153, 163, 845,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2330, 0, 3, 1922,
                                                                       671, 1958, 163, 173, 875,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2390, 0, 3, 1958,
                                                                       689, 1994, 173, 183, 905,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 2450, 0, 3, 2030,
                                                                       725, 2090, 203, 218, 935,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 2540, 0, 3, 2090,
                                                                       755, 2150, 218, 233, 980,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 2630, 0, 3, 2150,
                                                                       785, 2210, 233, 248, 1025,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 2720, 0, 3, 2210,
                                                                       815, 2270, 248, 263, 1070,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 2810, 0, 3, 2270,
                                                                       845, 2330, 263, 278, 1115,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 2900, 0, 3, 2330,
                                                                       875, 2390, 278, 293, 1160,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 2990, 0, 3, 2450,
                                                                       935, 2540, 323, 344, 1205,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 3116, 0, 3, 2540,
                                                                       980, 2630, 344, 365, 1268,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 3242, 0, 3, 2630,
                                                                       1025, 2720, 365, 386,
                                                                       1331, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 3368, 0, 3, 2720,
                                                                       1070, 2810, 386, 407,
                                                                       1394, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 3494, 0, 3, 2810,
                                                                       1115, 2900, 407, 428,
                                                                       1457, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3620, 3, 470, 473,
                                                                       1532, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3630, 3, 473, 476,
                                                                       1538, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3640, 3, 476, 479,
                                                                       1544, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3650, 3, 479, 482,
                                                                       1550, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3660, 3, 482, 485,
                                                                       1556, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3670, 3, 485, 488,
                                                                       1562, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3680, 3, 488, 491,
                                                                       1568, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3690, 3, 491, 494,
                                                                       1574, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 3700, 0, 3, 3620,
                                                                       1532, 3630, 500, 509,
                                                                       1616, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 3730, 0, 3, 3630,
                                                                       1538, 3640, 509, 518,
                                                                       1634, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 3760, 0, 3, 3640,
                                                                       1544, 3650, 518, 527,
                                                                       1652, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 3790, 0, 3, 3650,
                                                                       1550, 3660, 527, 536,
                                                                       1670, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 3820, 0, 3, 3660,
                                                                       1556, 3670, 536, 545,
                                                                       1688, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 3850, 0, 3, 3670,
                                                                       1562, 3680, 545, 554,
                                                                       1706, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 3880, 0, 3, 3680,
                                                                       1568, 3690, 554, 563,
                                                                       1724, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 3910, 0, 3, 3700,
                                                                       1616, 3730, 581, 599,
                                                                       1814, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 3970, 0, 3, 3730,
                                                                       1634, 3760, 599, 617,
                                                                       1850, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 4030, 0, 3, 3760,
                                                                       1652, 3790, 617, 635,
                                                                       1886, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 4090, 0, 3, 3790,
                                                                       1670, 3820, 635, 653,
                                                                       1922, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 4150, 0, 3, 3820,
                                                                       1688, 3850, 653, 671,
                                                                       1958, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 4210, 0, 3, 3850,
                                                                       1706, 3880, 671, 689,
                                                                       1994, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 4270, 0, 3, 3910,
                                                                       1814, 3970, 725, 755,
                                                                       2150, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 4370, 0, 3, 3970,
                                                                       1850, 4030, 755, 785,
                                                                       2210, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 4470, 0, 3, 4030,
                                                                       1886, 4090, 785, 815,
                                                                       2270, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 4570, 0, 3, 4090,
                                                                       1922, 4150, 815, 845,
                                                                       2330, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 4670, 0, 3, 4150,
                                                                       1958, 4210, 845, 875,
                                                                       2390, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 4770, 0, 3, 4270,
                                                                       2150, 4370, 935, 980,
                                                                       2630, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 4920, 0, 3, 4370,
                                                                       2210, 4470, 980, 1025,
                                                                       2720, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 5070, 0, 3, 4470,
                                                                       2270, 4570, 1025, 1070,
                                                                       2810, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 5220, 0, 3, 4570,
                                                                       2330, 4670, 1070, 1115,
                                                                       2900, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 5370, 0, 3, 4770,
                                                                       2630, 4920, 1205, 1268,
                                                                       3242, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 5580, 0, 3, 4920,
                                                                       2720, 5070, 1268, 1331,
                                                                       3368, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 5790, 0, 3, 5070,
                                                                       2810, 5220, 1331, 1394,
                                                                       3494, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 6000, 3, 1520,
                                                                       1526, 3620, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 6015, 3, 1526,
                                                                       1532, 3630, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 6030, 3, 1532,
                                                                       1538, 3640, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 6045, 3, 1538,
                                                                       1544, 3650, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 6060, 3, 1544,
                                                                       1550, 3660, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 6075, 3, 1550,
                                                                       1556, 3670, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 6090, 3, 1556,
                                                                       1562, 3680, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 6105, 3, 1562,
                                                                       1568, 3690, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 6120, 0, 3, 6000,
                                                                       3620, 6015, 1580, 1598,
                                                                       3700, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 6165, 0, 3, 6015,
                                                                       3630, 6030, 1598, 1616,
                                                                       3730, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 6210, 0, 3, 6030,
                                                                       3640, 6045, 1616, 1634,
                                                                       3760, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 6255, 0, 3, 6045,
                                                                       3650, 6060, 1634, 1652,
                                                                       3790, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 6300, 0, 3, 6060,
                                                                       3660, 6075, 1652, 1670,
                                                                       3820, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 6345, 0, 3, 6075,
                                                                       3670, 6090, 1670, 1688,
                                                                       3850, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 6390, 0, 3, 6090,
                                                                       3680, 6105, 1688, 1706,
                                                                       3880, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 6435, 0, 3, 6120,
                                                                       3700, 6165, 1742, 1778,
                                                                       3910, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 6525, 0, 3, 6165,
                                                                       3730, 6210, 1778, 1814,
                                                                       3970, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 6615, 0, 3, 6210,
                                                                       3760, 6255, 1814, 1850,
                                                                       4030, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 6705, 0, 3, 6255,
                                                                       3790, 6300, 1850, 1886,
                                                                       4090, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 6795, 0, 3, 6300,
                                                                       3820, 6345, 1886, 1922,
                                                                       4150, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 6885, 0, 3, 6345,
                                                                       3850, 6390, 1922, 1958,
                                                                       4210, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 6975, 0, 3, 6435,
                                                                       3910, 6525, 2030, 2090,
                                                                       4270, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 7125, 0, 3, 6525,
                                                                       3970, 6615, 2090, 2150,
                                                                       4370, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 7275, 0, 3, 6615,
                                                                       4030, 6705, 2150, 2210,
                                                                       4470, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 7425, 0, 3, 6705,
                                                                       4090, 6795, 2210, 2270,
                                                                       4570, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 7575, 0, 3, 6795,
                                                                       4150, 6885, 2270, 2330,
                                                                       4670, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 7725, 0, 3, 6975,
                                                                       4270, 7125, 2450, 2540,
                                                                       4770, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 7950, 0, 3, 7125,
                                                                       4370, 7275, 2540, 2630,
                                                                       4920, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 8175, 0, 3, 7275,
                                                                       4470, 7425, 2630, 2720,
                                                                       5070, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 8400, 0, 3, 7425,
                                                                       4570, 7575, 2720, 2810,
                                                                       5220, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 8625, 0, 3, 7725,
                                                                       4770, 7950, 2990, 3116,
                                                                       5370, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 8940, 0, 3, 7950,
                                                                       4920, 8175, 3116, 3242,
                                                                       5580, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 9255, 0, 3, 8175,
                                                                       5070, 8400, 3242, 3368,
                                                                       5790, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 9570, 3, 3620,
                                                                       3630, 6030, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 9591, 3, 3630,
                                                                       3640, 6045, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 9612, 3, 3640,
                                                                       3650, 6060, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 9633, 3, 3650,
                                                                       3660, 6075, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 9654, 3, 3660,
                                                                       3670, 6090, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 9675, 3, 3670,
                                                                       3680, 6105, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 9696, 0, 3, 9570,
                                                                       6030, 9591, 3700, 3730,
                                                                       6210, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 9759, 0, 3, 9591,
                                                                       6045, 9612, 3730, 3760,
                                                                       6255, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 9822, 0, 3, 9612,
                                                                       6060, 9633, 3760, 3790,
                                                                       6300, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 9885, 0, 3, 9633,
                                                                       6075, 9654, 3790, 3820,
                                                                       6345, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 9948, 0, 3, 9654,
                                                                       6090, 9675, 3820, 3850,
                                                                       6390, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 10011, 0, 3, 9696,
                                                                       6210, 9759, 3910, 3970,
                                                                       6615, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 10137, 0, 3, 9759,
                                                                       6255, 9822, 3970, 4030,
                                                                       6705, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 10263, 0, 3, 9822,
                                                                       6300, 9885, 4030, 4090,
                                                                       6795, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 10389, 0, 3, 9885,
                                                                       6345, 9948, 4090, 4150,
                                                                       6885, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 10515, 0, 3,
                                                                       10011, 6615, 10137, 4270,
                                                                       4370, 7275, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 10725, 0, 3,
                                                                       10137, 6705, 10263, 4370,
                                                                       4470, 7425, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 10935, 0, 3,
                                                                       10263, 6795, 10389, 4470,
                                                                       4570, 7575, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 11145, 0, 3,
                                                                       10515, 7275, 10725, 4770,
                                                                       4920, 8175, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 11460, 0, 3,
                                                                       10725, 7425, 10935, 4920,
                                                                       5070, 8400, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 11775, 0, 3,
                                                                       11145, 8175, 11460, 5370,
                                                                       5580, 9255, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 12216, 3, 6000,
                                                                       6015, 9570, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 12244, 3, 6015,
                                                                       6030, 9591, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 12272, 3, 6030,
                                                                       6045, 9612, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 12300, 3, 6045,
                                                                       6060, 9633, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 12328, 3, 6060,
                                                                       6075, 9654, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 12356, 3, 6075,
                                                                       6090, 9675, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 12384, 0, 3,
                                                                       12216, 9570, 12244, 6120,
                                                                       6165, 9696, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 12468, 0, 3,
                                                                       12244, 9591, 12272, 6165,
                                                                       6210, 9759, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 12552, 0, 3,
                                                                       12272, 9612, 12300, 6210,
                                                                       6255, 9822, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 12636, 0, 3,
                                                                       12300, 9633, 12328, 6255,
                                                                       6300, 9885, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 12720, 0, 3,
                                                                       12328, 9654, 12356, 6300,
                                                                       6345, 9948, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 12804, 0, 3,
                                                                       12384, 9696, 12468, 6435,
                                                                       6525, 10011, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 12972, 0, 3,
                                                                       12468, 9759, 12552, 6525,
                                                                       6615, 10137, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 13140, 0, 3,
                                                                       12552, 9822, 12636, 6615,
                                                                       6705, 10263, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 13308, 0, 3,
                                                                       12636, 9885, 12720, 6705,
                                                                       6795, 10389, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 13476, 0, 3,
                                                                       12804, 10011, 12972, 6975,
                                                                       7125, 10515, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 13756, 0, 3,
                                                                       12972, 10137, 13140, 7125,
                                                                       7275, 10725, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 14036, 0, 3,
                                                                       13140, 10263, 13308, 7275,
                                                                       7425, 10935, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 14316, 0, 3,
                                                                       13476, 10515, 13756, 7725,
                                                                       7950, 11145, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 14736, 0, 3,
                                                                       13756, 10725, 14036, 7950,
                                                                       8175, 11460, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 15156, 0, 3,
                                                                       14316, 11145, 14736, 8625,
                                                                       8940, 11775, ncols, gamma,
                                                                       p, q);

                    simdfunc::contract_primitives(buffer, 15744, 14316, 420, ncols);

                    simdfunc::contract_primitives(buffer, 16359, 15156, 588, ncols);
                }
            }
        }

        simdtrf::transform_i_inner(buffer, 16164, 15744, 15, 1, nmax);

        simdtrf::transform_i_inner(buffer, 16947, 16359, 21, 1, nmax);

        simdtrf::compute_hrr_gp_out_of_first(buffer, coordinates, 17220, 16164, 16947, 13,
                                             nmax);

        simdtrf::transform_p_inner(buffer, 17805, 17220, 15, 13, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 17805, 39, nmax);
    }

    for (size_t m = 0; m < 351; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
