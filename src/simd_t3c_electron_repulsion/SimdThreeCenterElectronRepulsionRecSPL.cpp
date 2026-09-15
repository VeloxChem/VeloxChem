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


#include "SimdThreeCenterElectronRepulsionRecSPL.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransformL.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_spl_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_spl_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 2053, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 51 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 2053, 1867, 135, dimensions);

        for (size_t i = 0; i < nprim_a; i++)
        {
            for (size_t j = 0; j < nprim_b; j++)
            {
                const auto p = a_exps[i] + b_exps[j];

                const auto mu = a_exps[i] * b_exps[j] / p;

                const auto fovl = a_norms[i] * b_norms[j];

                const auto fb = a_exps[i] / p;

                const auto fc = b_exps[j] / p;

                simdfunc::compute_pb(buffer, coordinates, 0, nmax, fb);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 7, 3, 9, ncols,
                                                             fj, 6, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 18, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 21, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 24, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 27, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 30, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 33, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 36, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 39, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 42, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 45, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 48, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 51, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 54, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 57, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 60, 3, 17, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 63, 3, 10, 18,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 72, 3, 11, 21,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 81, 3, 12, 24,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 90, 3, 13, 27,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 99, 3, 14, 30,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 108, 3, 15, 33,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 117, 3, 16, 36,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 126, 3, 8, 9, 39,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 132, 3, 9, 10, 42,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 138, 3, 10, 11,
                                                                       45, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 144, 3, 11, 12,
                                                                       48, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 150, 3, 12, 13,
                                                                       51, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 156, 3, 13, 14,
                                                                       54, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 162, 3, 14, 15,
                                                                       57, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 168, 3, 15, 16,
                                                                       60, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 174, 0, 3, 126,
                                                                       39, 132, 63, ncols, gamma,
                                                                       p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 192, 0, 3, 132,
                                                                       42, 138, 72, ncols, gamma,
                                                                       p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 210, 0, 3, 138,
                                                                       45, 144, 81, ncols, gamma,
                                                                       p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 228, 0, 3, 144,
                                                                       48, 150, 90, ncols, gamma,
                                                                       p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 246, 0, 3, 150,
                                                                       51, 156, 99, ncols, gamma,
                                                                       p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 264, 0, 3, 156,
                                                                       54, 162, 108, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 282, 0, 3, 162,
                                                                       57, 168, 117, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 300, 3, 39, 42,
                                                                       138, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 310, 3, 42, 45,
                                                                       144, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 320, 3, 45, 48,
                                                                       150, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 330, 3, 48, 51,
                                                                       156, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 340, 3, 51, 54,
                                                                       162, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 350, 3, 54, 57,
                                                                       168, ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 360, 0, 3, 300,
                                                                       138, 310, 210, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 390, 0, 3, 310,
                                                                       144, 320, 228, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 420, 0, 3, 320,
                                                                       150, 330, 246, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 450, 0, 3, 330,
                                                                       156, 340, 264, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 480, 0, 3, 340,
                                                                       162, 350, 282, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 510, 3, 126, 132,
                                                                       300, ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 525, 3, 132, 138,
                                                                       310, ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 540, 3, 138, 144,
                                                                       320, ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 555, 3, 144, 150,
                                                                       330, ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 570, 3, 150, 156,
                                                                       340, ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 585, 3, 156, 162,
                                                                       350, ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 600, 0, 3, 510,
                                                                       300, 525, 174, 192, 360,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 645, 0, 3, 525,
                                                                       310, 540, 192, 210, 390,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 690, 0, 3, 540,
                                                                       320, 555, 210, 228, 420,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 735, 0, 3, 555,
                                                                       330, 570, 228, 246, 450,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 780, 0, 3, 570,
                                                                       340, 585, 246, 264, 480,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 825, 3, 300, 310,
                                                                       540, ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 846, 3, 310, 320,
                                                                       555, ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 867, 3, 320, 330,
                                                                       570, ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 888, 3, 330, 340,
                                                                       585, ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 909, 0, 3, 825,
                                                                       540, 846, 360, 390, 690,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 972, 0, 3, 846,
                                                                       555, 867, 390, 420, 735,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 1035, 0, 3, 867,
                                                                       570, 888, 420, 450, 780,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 1098, 3, 510, 525,
                                                                       825, ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 1126, 3, 525, 540,
                                                                       846, ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 1154, 3, 540, 555,
                                                                       867, ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 1182, 3, 555, 570,
                                                                       888, ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 1210, 0, 3, 1098,
                                                                       825, 1126, 600, 645, 909,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 1294, 0, 3, 1126,
                                                                       846, 1154, 645, 690, 972,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 1378, 0, 3, 1154,
                                                                       867, 1182, 690, 735, 1035,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 1462, 3, 825, 846,
                                                                       1154, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 1498, 3, 846, 867,
                                                                       1182, ncols, gamma, p,
                                                                       q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 1534, 0, 3, 1462,
                                                                       1154, 1498, 909, 972,
                                                                       1378, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 1642, 3, 1098,
                                                                       1126, 1462, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 1687, 3, 1126,
                                                                       1154, 1498, ncols, gamma,
                                                                       p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 1732, 0, 3, 1642,
                                                                       1462, 1687, 1210, 1294,
                                                                       1534, ncols, gamma, p,
                                                                       q);

                    simdfunc::contract_primitives(buffer, 1867, 1732, 135, ncols);
                }
            }
        }

        simdtrf::transform_l_inner(buffer, 2002, 1867, 3, 1, nmax);

        simdtrf::transform_p_outer(values + n * npairs, nvalues, buffer, 2002, 17, nmax);
    }

    for (size_t m = 0; m < 51; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
