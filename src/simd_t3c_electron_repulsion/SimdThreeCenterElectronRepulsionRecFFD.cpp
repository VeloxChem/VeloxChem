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


#include "SimdThreeCenterElectronRepulsionRecFFD.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferFD.hpp"
#include "SimdTransferFF.hpp"
#include "SimdTransferFP.hpp"
#include "SimdTransferGD.hpp"
#include "SimdTransferGP.hpp"
#include "SimdTransferHP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformF.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_ffd_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_ffd_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 5380, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 245 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 5380, 2276, 674, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 7, 3, 8, ncols,
                                                             fj, 6, fq);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 386, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 389, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 392, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 395, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 398, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 401, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 404, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 407, 3, 10, 23,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 416, 3, 11, 26,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 425, 3, 12, 29,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 434, 3, 13, 32,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 443, 3, 14, 35,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 452, 3, 15, 38,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 461, 3, 23, 53,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 479, 3, 26, 59,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 497, 3, 29, 65,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 515, 3, 32, 71,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 533, 3, 35, 77,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 551, 3, 53, 103,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 581, 3, 59, 113,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 611, 3, 65, 123,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 641, 3, 71, 133,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 671, 3, 103, 173,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 716, 3, 113, 188,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 761, 3, 123, 203,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 806, 3, 173, 260,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 869, 3, 188, 281,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 932, 3, 260, 358,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1016, 3, 8, 9,
                                                                       386, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1022, 3, 9, 10,
                                                                       389, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1028, 3, 10, 11,
                                                                       392, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1034, 3, 11, 12,
                                                                       395, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1040, 3, 12, 13,
                                                                       398, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1046, 3, 13, 14,
                                                                       401, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1052, 3, 14, 15,
                                                                       404, ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1058, 0, 3, 1016,
                                                                       386, 1022, 407, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1076, 0, 3, 1022,
                                                                       389, 1028, 416, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1094, 0, 3, 1028,
                                                                       392, 1034, 425, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1112, 0, 3, 1034,
                                                                       395, 1040, 434, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1130, 0, 3, 1040,
                                                                       398, 1046, 443, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1148, 0, 3, 1046,
                                                                       401, 1052, 452, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1166, 0, 3, 1058,
                                                                       407, 1076, 41, 47, 461,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1202, 0, 3, 1076,
                                                                       416, 1094, 47, 53, 479,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1238, 0, 3, 1094,
                                                                       425, 1112, 53, 59, 497,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1274, 0, 3, 1112,
                                                                       434, 1130, 59, 65, 515,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1310, 0, 3, 1130,
                                                                       443, 1148, 65, 71, 533,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 1346, 0, 3, 1166,
                                                                       461, 1202, 83, 93, 551,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 1406, 0, 3, 1202,
                                                                       479, 1238, 93, 103, 581,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 1466, 0, 3, 1238,
                                                                       497, 1274, 103, 113, 611,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 1526, 0, 3, 1274,
                                                                       515, 1310, 113, 123, 641,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 1586, 0, 3, 1346,
                                                                       551, 1406, 143, 158, 671,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 1676, 0, 3, 1406,
                                                                       581, 1466, 158, 173, 716,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 1766, 0, 3, 1466,
                                                                       611, 1526, 173, 188, 761,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 1856, 0, 3, 1586,
                                                                       671, 1676, 218, 239, 806,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 1982, 0, 3, 1676,
                                                                       716, 1766, 239, 260, 869,
                                                                       ncols, gamma, p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 2108, 0, 3, 1856,
                                                                       806, 1982, 302, 330, 932,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 2276, 1346, 60, ncols);

                    simdfunc::contract_primitives(buffer, 2386, 1586, 90, ncols);

                    simdfunc::contract_primitives(buffer, 2551, 1856, 126, ncols);

                    simdfunc::contract_primitives(buffer, 2782, 2108, 168, ncols);
                }
            }
        }

        simdtrf::transform_d_inner(buffer, 2336, 2276, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 2476, 2386, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 2677, 2551, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 2950, 2782, 28, 1, nmax);

        simdtrf::compute_hrr_fp_out_of_first(buffer, coordinates, 3090, 2336, 2476, 5, nmax);

        simdtrf::compute_hrr_gp_out_of_first(buffer, coordinates, 3240, 2476, 2677, 5, nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 3465, 2677, 2950, 5, nmax);

        simdtrf::compute_hrr_fd_out_of_first(buffer, coordinates, 3780, 3090, 3240, 5, nmax);

        simdtrf::compute_hrr_gd_out_of_first(buffer, coordinates, 4080, 3240, 3465, 5, nmax);

        simdtrf::compute_hrr_ff_out_of_first(buffer, coordinates, 4530, 3780, 4080, 5, nmax);

        simdtrf::transform_f_inner(buffer, 5030, 4530, 10, 5, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 5030, 35, nmax);
    }

    for (size_t m = 0; m < 245; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
