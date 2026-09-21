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


#include "SimdNuclearPotentialRecGG.hpp"

#include <algorithm>
#include <cstddef>
#include <string>
#include <vector>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

#include "SimdNuclearPotentialVrrRecSD.hpp"
#include "SimdNuclearPotentialVrrRecSF.hpp"
#include "SimdNuclearPotentialVrrRecSG.hpp"
#include "SimdNuclearPotentialVrrRecSH.hpp"
#include "SimdNuclearPotentialVrrRecSI.hpp"
#include "SimdNuclearPotentialVrrRecSK.hpp"
#include "SimdNuclearPotentialVrrRecSL.hpp"
#include "SimdNuclearPotentialVrrRecSP.hpp"
#include "SimdTransferDG.hpp"
#include "SimdTransferDH.hpp"
#include "SimdTransferDI.hpp"
#include "SimdTransferFG.hpp"
#include "SimdTransferFH.hpp"
#include "SimdTransferGG.hpp"
#include "SimdTransferPG.hpp"
#include "SimdTransferPH.hpp"
#include "SimdTransferPI.hpp"
#include "SimdTransferPK.hpp"
#include "SimdTransformG.hpp"

namespace simdnpot {  // simdnpot namespace

auto
compute_gg_nuclear_potential(double                    *values,
                             const size_t               nvalues,
                             const CBasisFunction      &bra,
                             const CBasisFunction      &ket,
                             const CSimdMatrix         &coordinates,
                             const std::vector<double> &charges,
                             const std::vector<double> &points,
                             CSimdMatrix               &buffer,
                             const double               threshold) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_gg_nuclear_potential: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    errors::assertMsgCritical(
        points.size() == 3 * charges.size(),
        std::string("compute_gg_nuclear_potential: Expecting three coordinates for each charge"));

    if (charges.empty())
    {
        std::fill(values, values + 81 * nvalues, 0.0);

        return;
    }

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nprims = nprim_a * nprim_b;

    // NOTE: the pairs of primitives are screened with the threshold of the
    // integrals divided by their number and by the number of charges, as every
    // integral is a sum over both and the error of a sum is bounded by the
    // number of its terms.

    const auto terms = static_cast<double>(nprims * charges.size());

    const auto dimensions = simdfunc::make_column_dimensions(
        bra, ket, nvalues, coordinates, screenfunc::two_center_nuclear_potential_primitive_bound, threshold / terms);

    const auto nmax = simdfunc::prepare_buffer(buffer, 2052, 503, 145, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 81 * nvalues, 0.0);

        return;
    }

    errors::assertMsgCritical(dimensions.size() == nprim_a * nprim_b,
                              std::string("Dimensions do not match the pairs of primitives"));

    for (size_t i = 0; i < nprim_a; i++)
    {
        for (size_t j = 0; j < nprim_b; j++)
        {
            const auto ncols = dimensions[i * nprim_b + j];

            if (ncols == 0) continue;

            const auto p = a_exps[i] + b_exps[j];

            const auto mu = a_exps[i] * b_exps[j] / p;

            const auto fnpot = 2.0 * mathconst::pi_value() / p * a_norms[i] * b_norms[j];

            const auto fb = a_exps[i] / p;

            const auto fc = b_exps[j] / p;

            simdfunc::compute_pb(buffer, coordinates, 0, ncols, fb);

            simdfunc::compute_pair_exponent(buffer, coordinates, 6, ncols, mu);

            for (size_t ic = 0; ic < charges.size(); ic++)
            {
                const auto fz = fnpot * charges[ic];

                simdfunc::compute_pc(buffer, coordinates, 3, points, ic, ncols, fc);

                simdfunc::compute_full_npot_boys_function(buffer, coordinates, 7, 3, 8, ncols,
                                                          fz, 6, p);

                compute_prim_sp_nuclear_potential_0(buffer, 17, 0, 3, 8, 9, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 20, 0, 3, 9, 10, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 23, 0, 3, 10, 11, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 26, 0, 3, 11, 12, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 29, 0, 3, 12, 13, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 32, 0, 3, 13, 14, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 35, 0, 3, 14, 15, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 38, 0, 3, 15, 16, ncols);

                compute_prim_sd_nuclear_potential_0(buffer, 41, 0, 3, 8, 17, 9, 20, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 47, 0, 3, 9, 20, 10, 23, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 53, 0, 3, 10, 23, 11, 26, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 59, 0, 3, 11, 26, 12, 29, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 65, 0, 3, 12, 29, 13, 32, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 71, 0, 3, 13, 32, 14, 35, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 77, 0, 3, 14, 35, 15, 38, ncols, p);

                compute_prim_sf_nuclear_potential_0(buffer, 83, 0, 3, 17, 41, 20, 47, ncols, p);

                compute_prim_sf_nuclear_potential_0(buffer, 93, 0, 3, 20, 47, 23, 53, ncols, p);

                compute_prim_sf_nuclear_potential_0(buffer, 103, 0, 3, 23, 53, 26, 59, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 113, 0, 3, 26, 59, 29, 65, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 123, 0, 3, 29, 65, 32, 71, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 133, 0, 3, 32, 71, 35, 77, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 143, 0, 3, 41, 83, 47, 93, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 158, 0, 3, 47, 93, 53, 103, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 173, 0, 3, 53, 103, 59, 113, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 188, 0, 3, 59, 113, 65, 123, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 203, 0, 3, 65, 123, 71, 133, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 218, 0, 3, 83, 143, 93, 158, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 239, 0, 3, 93, 158, 103, 173, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 260, 0, 3, 103, 173, 113, 188, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 281, 0, 3, 113, 188, 123, 203, ncols,
                                                    p);

                compute_prim_si_nuclear_potential_0(buffer, 302, 0, 3, 143, 218, 158, 239, ncols,
                                                    p);

                compute_prim_si_nuclear_potential_0(buffer, 330, 0, 3, 158, 239, 173, 260, ncols,
                                                    p);

                compute_prim_si_nuclear_potential_0(buffer, 358, 0, 3, 173, 260, 188, 281, ncols,
                                                    p);

                compute_prim_sk_nuclear_potential_0(buffer, 386, 0, 3, 218, 302, 239, 330, ncols,
                                                    p);

                compute_prim_sk_nuclear_potential_0(buffer, 422, 0, 3, 239, 330, 260, 358, ncols,
                                                    p);

                compute_prim_sl_nuclear_potential_0(buffer, 458, 0, 3, 302, 386, 330, 422, ncols,
                                                    p);

                simdfunc::contract_primitives(buffer, 503, 143, 15, ncols);

                simdfunc::contract_primitives(buffer, 518, 218, 21, ncols);

                simdfunc::contract_primitives(buffer, 539, 302, 28, ncols);

                simdfunc::contract_primitives(buffer, 567, 386, 36, ncols);

                simdfunc::contract_primitives(buffer, 603, 458, 45, ncols);
            }
        }
    }

    simdtrf::compute_hrr_pg(buffer, coordinates, 648, 503, 518, 1, nmax);

    simdtrf::compute_hrr_ph(buffer, coordinates, 693, 518, 539, 1, nmax);

    simdtrf::compute_hrr_pi(buffer, coordinates, 756, 539, 567, 1, nmax);

    simdtrf::compute_hrr_pk(buffer, coordinates, 840, 567, 603, 1, nmax);

    simdtrf::compute_hrr_dg(buffer, coordinates, 948, 648, 693, 1, nmax);

    simdtrf::compute_hrr_dh(buffer, coordinates, 1038, 693, 756, 1, nmax);

    simdtrf::compute_hrr_di(buffer, coordinates, 1164, 756, 840, 1, nmax);

    simdtrf::compute_hrr_fg(buffer, coordinates, 1332, 948, 1038, 1, nmax);

    simdtrf::compute_hrr_fh(buffer, coordinates, 1482, 1038, 1164, 1, nmax);

    simdtrf::compute_hrr_gg(buffer, coordinates, 1692, 1332, 1482, 1, nmax);

    simdtrf::transform_g_inner(buffer, 1917, 1692, 15, 1, nmax);

    simdtrf::transform_g_outer(values, nvalues, buffer, 1917, 9, nmax);

    for (size_t m = 0; m < 81; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdnpot
