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


#include "SimdThreeCenterElectronRepulsionRecDPG.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferDP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_dpg_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_dpg_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 2381, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 135 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 2381, 1673, 294, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 7, 3, 7, ncols,
                                                             fj, 6, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 16, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 19, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 22, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 25, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 28, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 31, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 34, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 37, 0, 3, 8, 9,
                                                                       16, 19, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 43, 0, 3, 9, 10,
                                                                       19, 22, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 49, 0, 3, 10, 11,
                                                                       22, 25, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 55, 0, 3, 11, 12,
                                                                       25, 28, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 61, 0, 3, 12, 13,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 67, 0, 3, 13, 14,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 73, 0, 3, 16, 19,
                                                                       37, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 83, 0, 3, 19, 22,
                                                                       43, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 93, 0, 3, 22, 25,
                                                                       49, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 103, 0, 3, 25, 28,
                                                                       55, 61, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 113, 0, 3, 28, 31,
                                                                       61, 67, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 123, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 126, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 129, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 132, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 135, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 138, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 141, 3, 10, 22,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 150, 3, 11, 25,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 159, 3, 12, 28,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 168, 3, 13, 31,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 177, 3, 14, 34,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 186, 3, 22, 49,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 204, 3, 25, 55,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 222, 3, 28, 61,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 240, 3, 31, 67,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 258, 3, 49, 93,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 288, 3, 55, 103,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 318, 3, 61, 113,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 348, 3, 8, 9, 123,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 354, 3, 9, 10,
                                                                       126, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 360, 3, 10, 11,
                                                                       129, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 366, 3, 11, 12,
                                                                       132, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 372, 3, 12, 13,
                                                                       135, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 378, 3, 13, 14,
                                                                       138, ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 384, 0, 3, 348,
                                                                       123, 354, 141, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 402, 0, 3, 354,
                                                                       126, 360, 150, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 420, 0, 3, 360,
                                                                       129, 366, 159, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 438, 0, 3, 366,
                                                                       132, 372, 168, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 456, 0, 3, 372,
                                                                       135, 378, 177, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 474, 0, 3, 384,
                                                                       141, 402, 37, 43, 186,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 510, 0, 3, 402,
                                                                       150, 420, 43, 49, 204,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 546, 0, 3, 420,
                                                                       159, 438, 49, 55, 222,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 582, 0, 3, 438,
                                                                       168, 456, 55, 61, 240,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 618, 0, 3, 474,
                                                                       186, 510, 73, 83, 258,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 678, 0, 3, 510,
                                                                       204, 546, 83, 93, 288,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 738, 0, 3, 546,
                                                                       222, 582, 93, 103, 318,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 798, 3, 123, 126,
                                                                       360, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 808, 3, 126, 129,
                                                                       366, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 818, 3, 129, 132,
                                                                       372, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 828, 3, 132, 135,
                                                                       378, ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 838, 0, 3, 798,
                                                                       360, 808, 141, 150, 420,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 868, 0, 3, 808,
                                                                       366, 818, 150, 159, 438,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 898, 0, 3, 818,
                                                                       372, 828, 159, 168, 456,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 928, 0, 3, 838,
                                                                       420, 868, 186, 204, 546,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 988, 0, 3, 868,
                                                                       438, 898, 204, 222, 582,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 1048, 0, 3, 928,
                                                                       546, 988, 258, 288, 738,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 1148, 3, 348, 354,
                                                                       798, ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 1163, 3, 354, 360,
                                                                       808, ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 1178, 3, 360, 366,
                                                                       818, ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 1193, 3, 366, 372,
                                                                       828, ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 1208, 0, 3, 1148,
                                                                       798, 1163, 384, 402, 838,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 1253, 0, 3, 1163,
                                                                       808, 1178, 402, 420, 868,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 1298, 0, 3, 1178,
                                                                       818, 1193, 420, 438, 898,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 1343, 0, 3, 1208,
                                                                       838, 1253, 474, 510, 928,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 1433, 0, 3, 1253,
                                                                       868, 1298, 510, 546, 988,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 1523, 0, 3, 1343,
                                                                       928, 1433, 618, 678, 1048,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 1673, 1343, 90, ncols);

                    simdfunc::contract_primitives(buffer, 1817, 1523, 150, ncols);
                }
            }
        }

        simdtrf::transform_g_inner(buffer, 1763, 1673, 6, 1, nmax);

        simdtrf::transform_g_inner(buffer, 1967, 1817, 10, 1, nmax);

        simdtrf::compute_hrr_dp_out_of_first(buffer, coordinates, 2057, 1763, 1967, 9, nmax);

        simdtrf::transform_p_inner(buffer, 2219, 2057, 6, 9, nmax);

        simdtrf::transform_d_outer(values + n * npairs, nvalues, buffer, 2219, 27, nmax);
    }

    for (size_t m = 0; m < 135; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
