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


#include "SimdThreeCenterElectronRepulsionRecISK.hpp"

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
#include "SimdTransformI.hpp"
#include "SimdTransformK.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_isk_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_isk_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 40075, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 195 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 40075, 38647, 1008, dimensions);

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
                                                        5, 6, 7, 8, 9, 10, 11, 12, 13}, ncols,
                                                        fj, mu, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 20, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 23, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 26, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 29, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 32, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 35, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 38, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 41, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 44, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 47, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 50, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 53, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 56, 0, 3, 7, 8,
                                                                       20, 23, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 62, 0, 3, 8, 9,
                                                                       23, 26, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 68, 0, 3, 9, 10,
                                                                       26, 29, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 74, 0, 3, 10, 11,
                                                                       29, 32, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 80, 0, 3, 11, 12,
                                                                       32, 35, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 86, 0, 3, 12, 13,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 92, 0, 3, 13, 14,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 98, 0, 3, 14, 15,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 104, 0, 3, 15, 16,
                                                                       44, 47, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 110, 0, 3, 16, 17,
                                                                       47, 50, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 116, 0, 3, 17, 18,
                                                                       50, 53, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 122, 0, 3, 20, 23,
                                                                       56, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 132, 0, 3, 23, 26,
                                                                       62, 68, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 142, 0, 3, 26, 29,
                                                                       68, 74, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 152, 0, 3, 29, 32,
                                                                       74, 80, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 162, 0, 3, 32, 35,
                                                                       80, 86, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 172, 0, 3, 35, 38,
                                                                       86, 92, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 182, 0, 3, 38, 41,
                                                                       92, 98, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 192, 0, 3, 41, 44,
                                                                       98, 104, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 202, 0, 3, 44, 47,
                                                                       104, 110, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 212, 0, 3, 47, 50,
                                                                       110, 116, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 222, 0, 3, 56, 62,
                                                                       122, 132, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 237, 0, 3, 62, 68,
                                                                       132, 142, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 252, 0, 3, 68, 74,
                                                                       142, 152, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 267, 0, 3, 74, 80,
                                                                       152, 162, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 282, 0, 3, 80, 86,
                                                                       162, 172, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 297, 0, 3, 86, 92,
                                                                       172, 182, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 312, 0, 3, 92, 98,
                                                                       182, 192, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 327, 0, 3, 98,
                                                                       104, 192, 202, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 342, 0, 3, 104,
                                                                       110, 202, 212, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 357, 0, 3, 122,
                                                                       132, 222, 237, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 378, 0, 3, 132,
                                                                       142, 237, 252, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 399, 0, 3, 142,
                                                                       152, 252, 267, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 420, 0, 3, 152,
                                                                       162, 267, 282, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 441, 0, 3, 162,
                                                                       172, 282, 297, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 462, 0, 3, 172,
                                                                       182, 297, 312, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 483, 0, 3, 182,
                                                                       192, 312, 327, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 504, 0, 3, 192,
                                                                       202, 327, 342, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 525, 0, 3, 222,
                                                                       237, 357, 378, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 553, 0, 3, 237,
                                                                       252, 378, 399, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 581, 0, 3, 252,
                                                                       267, 399, 420, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 609, 0, 3, 267,
                                                                       282, 420, 441, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 637, 0, 3, 282,
                                                                       297, 441, 462, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 665, 0, 3, 297,
                                                                       312, 462, 483, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 693, 0, 3, 312,
                                                                       327, 483, 504, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 721, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 724, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 727, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 730, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 733, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 736, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 739, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 742, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 745, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 748, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 751, 3, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 754, 3, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 757, 3, 19, ncols,
                                                                       p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 760, 3, 7, 20,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 769, 3, 8, 23,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 778, 3, 9, 26,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 787, 3, 10, 29,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 796, 3, 11, 32,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 805, 3, 12, 35,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 814, 3, 13, 38,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 823, 3, 14, 41,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 832, 3, 15, 44,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 841, 3, 16, 47,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 850, 3, 17, 50,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 859, 3, 18, 53,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 868, 3, 20, 56,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 886, 3, 23, 62,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 904, 3, 26, 68,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 922, 3, 29, 74,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 940, 3, 32, 80,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 958, 3, 35, 86,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 976, 3, 38, 92,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 994, 3, 41, 98,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1012, 3, 44, 104,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1030, 3, 47, 110,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1048, 3, 50, 116,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1066, 3, 56, 122,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1096, 3, 62, 132,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1126, 3, 68, 142,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1156, 3, 74, 152,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1186, 3, 80, 162,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1216, 3, 86, 172,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1246, 3, 92, 182,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1276, 3, 98, 192,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1306, 3, 104, 202,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1336, 3, 110, 212,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1366, 3, 122, 222,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1411, 3, 132, 237,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1456, 3, 142, 252,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1501, 3, 152, 267,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1546, 3, 162, 282,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1591, 3, 172, 297,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1636, 3, 182, 312,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1681, 3, 192, 327,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1726, 3, 202, 342,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1771, 3, 222, 357,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1834, 3, 237, 378,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1897, 3, 252, 399,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1960, 3, 267, 420,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2023, 3, 282, 441,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2086, 3, 297, 462,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2149, 3, 312, 483,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2212, 3, 327, 504,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2275, 3, 357, 525,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2359, 3, 378, 553,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2443, 3, 399, 581,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2527, 3, 420, 609,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2611, 3, 441, 637,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2695, 3, 462, 665,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2779, 3, 483, 693,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2863, 3, 7, 8,
                                                                       727, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2869, 3, 8, 9,
                                                                       730, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2875, 3, 9, 10,
                                                                       733, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2881, 3, 10, 11,
                                                                       736, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2887, 3, 11, 12,
                                                                       739, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2893, 3, 12, 13,
                                                                       742, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2899, 3, 13, 14,
                                                                       745, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2905, 3, 14, 15,
                                                                       748, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2911, 3, 15, 16,
                                                                       751, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2917, 3, 16, 17,
                                                                       754, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2923, 3, 17, 18,
                                                                       757, ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2929, 0, 3, 2863,
                                                                       727, 2869, 778, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2947, 0, 3, 2869,
                                                                       730, 2875, 787, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2965, 0, 3, 2875,
                                                                       733, 2881, 796, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2983, 0, 3, 2881,
                                                                       736, 2887, 805, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3001, 0, 3, 2887,
                                                                       739, 2893, 814, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3019, 0, 3, 2893,
                                                                       742, 2899, 823, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3037, 0, 3, 2899,
                                                                       745, 2905, 832, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3055, 0, 3, 2905,
                                                                       748, 2911, 841, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3073, 0, 3, 2911,
                                                                       751, 2917, 850, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3091, 0, 3, 2917,
                                                                       754, 2923, 859, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3109, 0, 3, 2929,
                                                                       778, 2947, 56, 62, 904,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3145, 0, 3, 2947,
                                                                       787, 2965, 62, 68, 922,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3181, 0, 3, 2965,
                                                                       796, 2983, 68, 74, 940,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3217, 0, 3, 2983,
                                                                       805, 3001, 74, 80, 958,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3253, 0, 3, 3001,
                                                                       814, 3019, 80, 86, 976,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3289, 0, 3, 3019,
                                                                       823, 3037, 86, 92, 994,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3325, 0, 3, 3037,
                                                                       832, 3055, 92, 98, 1012,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3361, 0, 3, 3055,
                                                                       841, 3073, 98, 104, 1030,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3397, 0, 3, 3073,
                                                                       850, 3091, 104, 110, 1048,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3433, 0, 3, 3109,
                                                                       904, 3145, 122, 132, 1126,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3493, 0, 3, 3145,
                                                                       922, 3181, 132, 142, 1156,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3553, 0, 3, 3181,
                                                                       940, 3217, 142, 152, 1186,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3613, 0, 3, 3217,
                                                                       958, 3253, 152, 162, 1216,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3673, 0, 3, 3253,
                                                                       976, 3289, 162, 172, 1246,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3733, 0, 3, 3289,
                                                                       994, 3325, 172, 182, 1276,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3793, 0, 3, 3325,
                                                                       1012, 3361, 182, 192,
                                                                       1306, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3853, 0, 3, 3361,
                                                                       1030, 3397, 192, 202,
                                                                       1336, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 3913, 0, 3, 3433,
                                                                       1126, 3493, 222, 237,
                                                                       1456, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4003, 0, 3, 3493,
                                                                       1156, 3553, 237, 252,
                                                                       1501, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4093, 0, 3, 3553,
                                                                       1186, 3613, 252, 267,
                                                                       1546, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4183, 0, 3, 3613,
                                                                       1216, 3673, 267, 282,
                                                                       1591, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4273, 0, 3, 3673,
                                                                       1246, 3733, 282, 297,
                                                                       1636, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4363, 0, 3, 3733,
                                                                       1276, 3793, 297, 312,
                                                                       1681, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4453, 0, 3, 3793,
                                                                       1306, 3853, 312, 327,
                                                                       1726, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 4543, 0, 3, 3913,
                                                                       1456, 4003, 357, 378,
                                                                       1897, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 4669, 0, 3, 4003,
                                                                       1501, 4093, 378, 399,
                                                                       1960, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 4795, 0, 3, 4093,
                                                                       1546, 4183, 399, 420,
                                                                       2023, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 4921, 0, 3, 4183,
                                                                       1591, 4273, 420, 441,
                                                                       2086, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 5047, 0, 3, 4273,
                                                                       1636, 4363, 441, 462,
                                                                       2149, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 5173, 0, 3, 4363,
                                                                       1681, 4453, 462, 483,
                                                                       2212, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 5299, 0, 3, 4543,
                                                                       1897, 4669, 525, 553,
                                                                       2443, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 5467, 0, 3, 4669,
                                                                       1960, 4795, 553, 581,
                                                                       2527, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 5635, 0, 3, 4795,
                                                                       2023, 4921, 581, 609,
                                                                       2611, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 5803, 0, 3, 4921,
                                                                       2086, 5047, 609, 637,
                                                                       2695, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 5971, 0, 3, 5047,
                                                                       2149, 5173, 637, 665,
                                                                       2779, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6139, 3, 721, 724,
                                                                       2863, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6149, 3, 724, 727,
                                                                       2869, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6159, 3, 727, 730,
                                                                       2875, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6169, 3, 730, 733,
                                                                       2881, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6179, 3, 733, 736,
                                                                       2887, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6189, 3, 736, 739,
                                                                       2893, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6199, 3, 739, 742,
                                                                       2899, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6209, 3, 742, 745,
                                                                       2905, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6219, 3, 745, 748,
                                                                       2911, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6229, 3, 748, 751,
                                                                       2917, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6239, 3, 751, 754,
                                                                       2923, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6249, 0, 3, 6139,
                                                                       2863, 6149, 760, 769,
                                                                       2929, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6279, 0, 3, 6149,
                                                                       2869, 6159, 769, 778,
                                                                       2947, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6309, 0, 3, 6159,
                                                                       2875, 6169, 778, 787,
                                                                       2965, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6339, 0, 3, 6169,
                                                                       2881, 6179, 787, 796,
                                                                       2983, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6369, 0, 3, 6179,
                                                                       2887, 6189, 796, 805,
                                                                       3001, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6399, 0, 3, 6189,
                                                                       2893, 6199, 805, 814,
                                                                       3019, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6429, 0, 3, 6199,
                                                                       2899, 6209, 814, 823,
                                                                       3037, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6459, 0, 3, 6209,
                                                                       2905, 6219, 823, 832,
                                                                       3055, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6489, 0, 3, 6219,
                                                                       2911, 6229, 832, 841,
                                                                       3073, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6519, 0, 3, 6229,
                                                                       2917, 6239, 841, 850,
                                                                       3091, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6549, 0, 3, 6249,
                                                                       2929, 6279, 868, 886,
                                                                       3109, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6609, 0, 3, 6279,
                                                                       2947, 6309, 886, 904,
                                                                       3145, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6669, 0, 3, 6309,
                                                                       2965, 6339, 904, 922,
                                                                       3181, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6729, 0, 3, 6339,
                                                                       2983, 6369, 922, 940,
                                                                       3217, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6789, 0, 3, 6369,
                                                                       3001, 6399, 940, 958,
                                                                       3253, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6849, 0, 3, 6399,
                                                                       3019, 6429, 958, 976,
                                                                       3289, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6909, 0, 3, 6429,
                                                                       3037, 6459, 976, 994,
                                                                       3325, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6969, 0, 3, 6459,
                                                                       3055, 6489, 994, 1012,
                                                                       3361, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 7029, 0, 3, 6489,
                                                                       3073, 6519, 1012, 1030,
                                                                       3397, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7089, 0, 3, 6549,
                                                                       3109, 6609, 1066, 1096,
                                                                       3433, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7189, 0, 3, 6609,
                                                                       3145, 6669, 1096, 1126,
                                                                       3493, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7289, 0, 3, 6669,
                                                                       3181, 6729, 1126, 1156,
                                                                       3553, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7389, 0, 3, 6729,
                                                                       3217, 6789, 1156, 1186,
                                                                       3613, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7489, 0, 3, 6789,
                                                                       3253, 6849, 1186, 1216,
                                                                       3673, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7589, 0, 3, 6849,
                                                                       3289, 6909, 1216, 1246,
                                                                       3733, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7689, 0, 3, 6909,
                                                                       3325, 6969, 1246, 1276,
                                                                       3793, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7789, 0, 3, 6969,
                                                                       3361, 7029, 1276, 1306,
                                                                       3853, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 7889, 0, 3, 7089,
                                                                       3433, 7189, 1366, 1411,
                                                                       3913, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 8039, 0, 3, 7189,
                                                                       3493, 7289, 1411, 1456,
                                                                       4003, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 8189, 0, 3, 7289,
                                                                       3553, 7389, 1456, 1501,
                                                                       4093, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 8339, 0, 3, 7389,
                                                                       3613, 7489, 1501, 1546,
                                                                       4183, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 8489, 0, 3, 7489,
                                                                       3673, 7589, 1546, 1591,
                                                                       4273, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 8639, 0, 3, 7589,
                                                                       3733, 7689, 1591, 1636,
                                                                       4363, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 8789, 0, 3, 7689,
                                                                       3793, 7789, 1636, 1681,
                                                                       4453, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 8939, 0, 3, 7889,
                                                                       3913, 8039, 1771, 1834,
                                                                       4543, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 9149, 0, 3, 8039,
                                                                       4003, 8189, 1834, 1897,
                                                                       4669, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 9359, 0, 3, 8189,
                                                                       4093, 8339, 1897, 1960,
                                                                       4795, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 9569, 0, 3, 8339,
                                                                       4183, 8489, 1960, 2023,
                                                                       4921, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 9779, 0, 3, 8489,
                                                                       4273, 8639, 2023, 2086,
                                                                       5047, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 9989, 0, 3, 8639,
                                                                       4363, 8789, 2086, 2149,
                                                                       5173, ncols, gamma, p,
                                                                       q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 10199, 0, 3, 8939,
                                                                       4543, 9149, 2275, 2359,
                                                                       5299, ncols, gamma, p,
                                                                       q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 10479, 0, 3, 9149,
                                                                       4669, 9359, 2359, 2443,
                                                                       5467, ncols, gamma, p,
                                                                       q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 10759, 0, 3, 9359,
                                                                       4795, 9569, 2443, 2527,
                                                                       5635, ncols, gamma, p,
                                                                       q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 11039, 0, 3, 9569,
                                                                       4921, 9779, 2527, 2611,
                                                                       5803, ncols, gamma, p,
                                                                       q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 11319, 0, 3, 9779,
                                                                       5047, 9989, 2611, 2695,
                                                                       5971, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 11599, 3, 2863,
                                                                       2869, 6159, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 11614, 3, 2869,
                                                                       2875, 6169, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 11629, 3, 2875,
                                                                       2881, 6179, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 11644, 3, 2881,
                                                                       2887, 6189, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 11659, 3, 2887,
                                                                       2893, 6199, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 11674, 3, 2893,
                                                                       2899, 6209, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 11689, 3, 2899,
                                                                       2905, 6219, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 11704, 3, 2905,
                                                                       2911, 6229, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 11719, 3, 2911,
                                                                       2917, 6239, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 11734, 0, 3,
                                                                       11599, 6159, 11614, 2929,
                                                                       2947, 6309, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 11779, 0, 3,
                                                                       11614, 6169, 11629, 2947,
                                                                       2965, 6339, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 11824, 0, 3,
                                                                       11629, 6179, 11644, 2965,
                                                                       2983, 6369, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 11869, 0, 3,
                                                                       11644, 6189, 11659, 2983,
                                                                       3001, 6399, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 11914, 0, 3,
                                                                       11659, 6199, 11674, 3001,
                                                                       3019, 6429, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 11959, 0, 3,
                                                                       11674, 6209, 11689, 3019,
                                                                       3037, 6459, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 12004, 0, 3,
                                                                       11689, 6219, 11704, 3037,
                                                                       3055, 6489, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 12049, 0, 3,
                                                                       11704, 6229, 11719, 3055,
                                                                       3073, 6519, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 12094, 0, 3,
                                                                       11734, 6309, 11779, 3109,
                                                                       3145, 6669, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 12184, 0, 3,
                                                                       11779, 6339, 11824, 3145,
                                                                       3181, 6729, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 12274, 0, 3,
                                                                       11824, 6369, 11869, 3181,
                                                                       3217, 6789, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 12364, 0, 3,
                                                                       11869, 6399, 11914, 3217,
                                                                       3253, 6849, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 12454, 0, 3,
                                                                       11914, 6429, 11959, 3253,
                                                                       3289, 6909, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 12544, 0, 3,
                                                                       11959, 6459, 12004, 3289,
                                                                       3325, 6969, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 12634, 0, 3,
                                                                       12004, 6489, 12049, 3325,
                                                                       3361, 7029, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 12724, 0, 3,
                                                                       12094, 6669, 12184, 3433,
                                                                       3493, 7289, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 12874, 0, 3,
                                                                       12184, 6729, 12274, 3493,
                                                                       3553, 7389, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 13024, 0, 3,
                                                                       12274, 6789, 12364, 3553,
                                                                       3613, 7489, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 13174, 0, 3,
                                                                       12364, 6849, 12454, 3613,
                                                                       3673, 7589, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 13324, 0, 3,
                                                                       12454, 6909, 12544, 3673,
                                                                       3733, 7689, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 13474, 0, 3,
                                                                       12544, 6969, 12634, 3733,
                                                                       3793, 7789, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 13624, 0, 3,
                                                                       12724, 7289, 12874, 3913,
                                                                       4003, 8189, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 13849, 0, 3,
                                                                       12874, 7389, 13024, 4003,
                                                                       4093, 8339, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 14074, 0, 3,
                                                                       13024, 7489, 13174, 4093,
                                                                       4183, 8489, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 14299, 0, 3,
                                                                       13174, 7589, 13324, 4183,
                                                                       4273, 8639, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 14524, 0, 3,
                                                                       13324, 7689, 13474, 4273,
                                                                       4363, 8789, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 14749, 0, 3,
                                                                       13624, 8189, 13849, 4543,
                                                                       4669, 9359, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 15064, 0, 3,
                                                                       13849, 8339, 14074, 4669,
                                                                       4795, 9569, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 15379, 0, 3,
                                                                       14074, 8489, 14299, 4795,
                                                                       4921, 9779, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 15694, 0, 3,
                                                                       14299, 8639, 14524, 4921,
                                                                       5047, 9989, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 16009, 0, 3,
                                                                       14749, 9359, 15064, 5299,
                                                                       5467, 10759, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 16429, 0, 3,
                                                                       15064, 9569, 15379, 5467,
                                                                       5635, 11039, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 16849, 0, 3,
                                                                       15379, 9779, 15694, 5635,
                                                                       5803, 11319, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17269, 3, 6139,
                                                                       6149, 11599, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17290, 3, 6149,
                                                                       6159, 11614, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17311, 3, 6159,
                                                                       6169, 11629, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17332, 3, 6169,
                                                                       6179, 11644, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17353, 3, 6179,
                                                                       6189, 11659, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17374, 3, 6189,
                                                                       6199, 11674, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17395, 3, 6199,
                                                                       6209, 11689, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17416, 3, 6209,
                                                                       6219, 11704, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17437, 3, 6219,
                                                                       6229, 11719, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 17458, 0, 3,
                                                                       17269, 11599, 17290, 6249,
                                                                       6279, 11734, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 17521, 0, 3,
                                                                       17290, 11614, 17311, 6279,
                                                                       6309, 11779, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 17584, 0, 3,
                                                                       17311, 11629, 17332, 6309,
                                                                       6339, 11824, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 17647, 0, 3,
                                                                       17332, 11644, 17353, 6339,
                                                                       6369, 11869, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 17710, 0, 3,
                                                                       17353, 11659, 17374, 6369,
                                                                       6399, 11914, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 17773, 0, 3,
                                                                       17374, 11674, 17395, 6399,
                                                                       6429, 11959, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 17836, 0, 3,
                                                                       17395, 11689, 17416, 6429,
                                                                       6459, 12004, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 17899, 0, 3,
                                                                       17416, 11704, 17437, 6459,
                                                                       6489, 12049, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 17962, 0, 3,
                                                                       17458, 11734, 17521, 6549,
                                                                       6609, 12094, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 18088, 0, 3,
                                                                       17521, 11779, 17584, 6609,
                                                                       6669, 12184, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 18214, 0, 3,
                                                                       17584, 11824, 17647, 6669,
                                                                       6729, 12274, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 18340, 0, 3,
                                                                       17647, 11869, 17710, 6729,
                                                                       6789, 12364, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 18466, 0, 3,
                                                                       17710, 11914, 17773, 6789,
                                                                       6849, 12454, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 18592, 0, 3,
                                                                       17773, 11959, 17836, 6849,
                                                                       6909, 12544, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 18718, 0, 3,
                                                                       17836, 12004, 17899, 6909,
                                                                       6969, 12634, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 18844, 0, 3,
                                                                       17962, 12094, 18088, 7089,
                                                                       7189, 12724, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 19054, 0, 3,
                                                                       18088, 12184, 18214, 7189,
                                                                       7289, 12874, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 19264, 0, 3,
                                                                       18214, 12274, 18340, 7289,
                                                                       7389, 13024, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 19474, 0, 3,
                                                                       18340, 12364, 18466, 7389,
                                                                       7489, 13174, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 19684, 0, 3,
                                                                       18466, 12454, 18592, 7489,
                                                                       7589, 13324, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 19894, 0, 3,
                                                                       18592, 12544, 18718, 7589,
                                                                       7689, 13474, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 20104, 0, 3,
                                                                       18844, 12724, 19054, 7889,
                                                                       8039, 13624, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 20419, 0, 3,
                                                                       19054, 12874, 19264, 8039,
                                                                       8189, 13849, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 20734, 0, 3,
                                                                       19264, 13024, 19474, 8189,
                                                                       8339, 14074, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 21049, 0, 3,
                                                                       19474, 13174, 19684, 8339,
                                                                       8489, 14299, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 21364, 0, 3,
                                                                       19684, 13324, 19894, 8489,
                                                                       8639, 14524, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 21679, 0, 3,
                                                                       20104, 13624, 20419, 8939,
                                                                       9149, 14749, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 22120, 0, 3,
                                                                       20419, 13849, 20734, 9149,
                                                                       9359, 15064, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 22561, 0, 3,
                                                                       20734, 14074, 21049, 9359,
                                                                       9569, 15379, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 23002, 0, 3,
                                                                       21049, 14299, 21364, 9569,
                                                                       9779, 15694, ncols, gamma,
                                                                       p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 23443, 0, 3,
                                                                       21679, 14749, 22120,
                                                                       10199, 10479, 16009,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 24031, 0, 3,
                                                                       22120, 15064, 22561,
                                                                       10479, 10759, 16429,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 24619, 0, 3,
                                                                       22561, 15379, 23002,
                                                                       10759, 11039, 16849,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 25207, 3, 11599,
                                                                       11614, 17311, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 25235, 3, 11614,
                                                                       11629, 17332, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 25263, 3, 11629,
                                                                       11644, 17353, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 25291, 3, 11644,
                                                                       11659, 17374, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 25319, 3, 11659,
                                                                       11674, 17395, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 25347, 3, 11674,
                                                                       11689, 17416, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 25375, 3, 11689,
                                                                       11704, 17437, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 25403, 0, 3,
                                                                       25207, 17311, 25235,
                                                                       11734, 11779, 17584,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 25487, 0, 3,
                                                                       25235, 17332, 25263,
                                                                       11779, 11824, 17647,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 25571, 0, 3,
                                                                       25263, 17353, 25291,
                                                                       11824, 11869, 17710,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 25655, 0, 3,
                                                                       25291, 17374, 25319,
                                                                       11869, 11914, 17773,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 25739, 0, 3,
                                                                       25319, 17395, 25347,
                                                                       11914, 11959, 17836,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 25823, 0, 3,
                                                                       25347, 17416, 25375,
                                                                       11959, 12004, 17899,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 25907, 0, 3,
                                                                       25403, 17584, 25487,
                                                                       12094, 12184, 18214,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 26075, 0, 3,
                                                                       25487, 17647, 25571,
                                                                       12184, 12274, 18340,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 26243, 0, 3,
                                                                       25571, 17710, 25655,
                                                                       12274, 12364, 18466,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 26411, 0, 3,
                                                                       25655, 17773, 25739,
                                                                       12364, 12454, 18592,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 26579, 0, 3,
                                                                       25739, 17836, 25823,
                                                                       12454, 12544, 18718,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 26747, 0, 3,
                                                                       25907, 18214, 26075,
                                                                       12724, 12874, 19264,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 27027, 0, 3,
                                                                       26075, 18340, 26243,
                                                                       12874, 13024, 19474,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 27307, 0, 3,
                                                                       26243, 18466, 26411,
                                                                       13024, 13174, 19684,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 27587, 0, 3,
                                                                       26411, 18592, 26579,
                                                                       13174, 13324, 19894,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 27867, 0, 3,
                                                                       26747, 19264, 27027,
                                                                       13624, 13849, 20734,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 28287, 0, 3,
                                                                       27027, 19474, 27307,
                                                                       13849, 14074, 21049,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 28707, 0, 3,
                                                                       27307, 19684, 27587,
                                                                       14074, 14299, 21364,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 29127, 0, 3,
                                                                       27867, 20734, 28287,
                                                                       14749, 15064, 22561,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 29715, 0, 3,
                                                                       28287, 21049, 28707,
                                                                       15064, 15379, 23002,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 30303, 0, 3,
                                                                       29127, 22561, 29715,
                                                                       16009, 16429, 24619,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 31087, 3, 17269,
                                                                       17290, 25207, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 31123, 3, 17290,
                                                                       17311, 25235, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 31159, 3, 17311,
                                                                       17332, 25263, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 31195, 3, 17332,
                                                                       17353, 25291, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 31231, 3, 17353,
                                                                       17374, 25319, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 31267, 3, 17374,
                                                                       17395, 25347, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 31303, 3, 17395,
                                                                       17416, 25375, ncols,
                                                                       gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 31339, 0, 3,
                                                                       31087, 25207, 31123,
                                                                       17458, 17521, 25403,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 31447, 0, 3,
                                                                       31123, 25235, 31159,
                                                                       17521, 17584, 25487,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 31555, 0, 3,
                                                                       31159, 25263, 31195,
                                                                       17584, 17647, 25571,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 31663, 0, 3,
                                                                       31195, 25291, 31231,
                                                                       17647, 17710, 25655,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 31771, 0, 3,
                                                                       31231, 25319, 31267,
                                                                       17710, 17773, 25739,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 31879, 0, 3,
                                                                       31267, 25347, 31303,
                                                                       17773, 17836, 25823,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 31987, 0, 3,
                                                                       31339, 25403, 31447,
                                                                       17962, 18088, 25907,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 32203, 0, 3,
                                                                       31447, 25487, 31555,
                                                                       18088, 18214, 26075,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 32419, 0, 3,
                                                                       31555, 25571, 31663,
                                                                       18214, 18340, 26243,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 32635, 0, 3,
                                                                       31663, 25655, 31771,
                                                                       18340, 18466, 26411,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 32851, 0, 3,
                                                                       31771, 25739, 31879,
                                                                       18466, 18592, 26579,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 33067, 0, 3,
                                                                       31987, 25907, 32203,
                                                                       18844, 19054, 26747,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 33427, 0, 3,
                                                                       32203, 26075, 32419,
                                                                       19054, 19264, 27027,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 33787, 0, 3,
                                                                       32419, 26243, 32635,
                                                                       19264, 19474, 27307,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 34147, 0, 3,
                                                                       32635, 26411, 32851,
                                                                       19474, 19684, 27587,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 34507, 0, 3,
                                                                       33067, 26747, 33427,
                                                                       20104, 20419, 27867,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 35047, 0, 3,
                                                                       33427, 27027, 33787,
                                                                       20419, 20734, 28287,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 35587, 0, 3,
                                                                       33787, 27307, 34147,
                                                                       20734, 21049, 28707,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 36127, 0, 3,
                                                                       34507, 27867, 35047,
                                                                       21679, 22120, 29127,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 36883, 0, 3,
                                                                       35047, 28287, 35587,
                                                                       22120, 22561, 29715,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 37639, 0, 3,
                                                                       36127, 29127, 36883,
                                                                       23443, 24031, 30303,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 38647, 37639, 1008, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 39655, 38647, 28, 1, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 39655, 15, nmax);
    }

    for (size_t m = 0; m < 195; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
