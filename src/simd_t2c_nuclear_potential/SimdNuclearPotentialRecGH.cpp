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


#include "SimdNuclearPotentialRecGH.hpp"

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
#include "SimdNuclearPotentialVrrRecSM.hpp"
#include "SimdNuclearPotentialVrrRecSP.hpp"
#include "SimdTransferDH.hpp"
#include "SimdTransferDI.hpp"
#include "SimdTransferDK.hpp"
#include "SimdTransferFH.hpp"
#include "SimdTransferFI.hpp"
#include "SimdTransferGH.hpp"
#include "SimdTransferPH.hpp"
#include "SimdTransferPI.hpp"
#include "SimdTransferPK.hpp"
#include "SimdTransferPL.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformH.hpp"

namespace simdnpot {  // simdnpot namespace

auto
compute_gh_nuclear_potential(double                    *values,
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
            false, std::string("compute_gh_nuclear_potential: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    errors::assertMsgCritical(
        points.size() == 3 * charges.size(),
        std::string("compute_gh_nuclear_potential: Expecting three coordinates for each charge"));

    if (charges.empty())
    {
        std::fill(values, values + 99 * nvalues, 0.0);

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

    const auto nmax = simdfunc::prepare_buffer(buffer, 2778, 723, 185, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 99 * nvalues, 0.0);

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

                simdfunc::compute_full_npot_boys_function(buffer, coordinates, 7, 3, 9, ncols,
                                                          fz, 6, p);

                compute_prim_sp_nuclear_potential_0(buffer, 18, 0, 3, 8, 9, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 21, 0, 3, 9, 10, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 24, 0, 3, 10, 11, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 27, 0, 3, 11, 12, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 30, 0, 3, 12, 13, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 33, 0, 3, 13, 14, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 36, 0, 3, 14, 15, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 39, 0, 3, 15, 16, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 42, 0, 3, 16, 17, ncols);

                compute_prim_sd_nuclear_potential_0(buffer, 45, 0, 3, 8, 18, 9, 21, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 51, 0, 3, 9, 21, 10, 24, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 57, 0, 3, 10, 24, 11, 27, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 63, 0, 3, 11, 27, 12, 30, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 69, 0, 3, 12, 30, 13, 33, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 75, 0, 3, 13, 33, 14, 36, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 81, 0, 3, 14, 36, 15, 39, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 87, 0, 3, 15, 39, 16, 42, ncols, p);

                compute_prim_sf_nuclear_potential_0(buffer, 93, 0, 3, 18, 45, 21, 51, ncols, p);

                compute_prim_sf_nuclear_potential_0(buffer, 103, 0, 3, 21, 51, 24, 57, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 113, 0, 3, 24, 57, 27, 63, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 123, 0, 3, 27, 63, 30, 69, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 133, 0, 3, 30, 69, 33, 75, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 143, 0, 3, 33, 75, 36, 81, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 153, 0, 3, 36, 81, 39, 87, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 163, 0, 3, 45, 93, 51, 103, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 178, 0, 3, 51, 103, 57, 113, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 193, 0, 3, 57, 113, 63, 123, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 208, 0, 3, 63, 123, 69, 133, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 223, 0, 3, 69, 133, 75, 143, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 238, 0, 3, 75, 143, 81, 153, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 253, 0, 3, 93, 163, 103, 178, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 274, 0, 3, 103, 178, 113, 193, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 295, 0, 3, 113, 193, 123, 208, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 316, 0, 3, 123, 208, 133, 223, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 337, 0, 3, 133, 223, 143, 238, ncols,
                                                    p);

                compute_prim_si_nuclear_potential_0(buffer, 358, 0, 3, 163, 253, 178, 274, ncols,
                                                    p);

                compute_prim_si_nuclear_potential_0(buffer, 386, 0, 3, 178, 274, 193, 295, ncols,
                                                    p);

                compute_prim_si_nuclear_potential_0(buffer, 414, 0, 3, 193, 295, 208, 316, ncols,
                                                    p);

                compute_prim_si_nuclear_potential_0(buffer, 442, 0, 3, 208, 316, 223, 337, ncols,
                                                    p);

                compute_prim_sk_nuclear_potential_0(buffer, 470, 0, 3, 253, 358, 274, 386, ncols,
                                                    p);

                compute_prim_sk_nuclear_potential_0(buffer, 506, 0, 3, 274, 386, 295, 414, ncols,
                                                    p);

                compute_prim_sk_nuclear_potential_0(buffer, 542, 0, 3, 295, 414, 316, 442, ncols,
                                                    p);

                compute_prim_sl_nuclear_potential_0(buffer, 578, 0, 3, 358, 470, 386, 506, ncols,
                                                    p);

                compute_prim_sl_nuclear_potential_0(buffer, 623, 0, 3, 386, 506, 414, 542, ncols,
                                                    p);

                compute_prim_sm_nuclear_potential_0(buffer, 668, 0, 3, 470, 578, 506, 623, ncols,
                                                    p);

                simdfunc::contract_primitives(buffer, 723, 253, 21, ncols);

                simdfunc::contract_primitives(buffer, 744, 358, 28, ncols);

                simdfunc::contract_primitives(buffer, 772, 470, 36, ncols);

                simdfunc::contract_primitives(buffer, 808, 578, 45, ncols);

                simdfunc::contract_primitives(buffer, 853, 668, 55, ncols);
            }
        }
    }

    simdtrf::compute_hrr_ph(buffer, coordinates, 908, 723, 744, 1, nmax);

    simdtrf::compute_hrr_pi(buffer, coordinates, 971, 744, 772, 1, nmax);

    simdtrf::compute_hrr_pk(buffer, coordinates, 1055, 772, 808, 1, nmax);

    simdtrf::compute_hrr_pl(buffer, coordinates, 1163, 808, 853, 1, nmax);

    simdtrf::compute_hrr_dh(buffer, coordinates, 1298, 908, 971, 1, nmax);

    simdtrf::compute_hrr_di(buffer, coordinates, 1424, 971, 1055, 1, nmax);

    simdtrf::compute_hrr_dk(buffer, coordinates, 1592, 1055, 1163, 1, nmax);

    simdtrf::compute_hrr_fh(buffer, coordinates, 1808, 1298, 1424, 1, nmax);

    simdtrf::compute_hrr_fi(buffer, coordinates, 2018, 1424, 1592, 1, nmax);

    simdtrf::compute_hrr_gh(buffer, coordinates, 2298, 1808, 2018, 1, nmax);

    simdtrf::transform_h_inner(buffer, 2613, 2298, 15, 1, nmax);

    simdtrf::transform_g_outer(values, nvalues, buffer, 2613, 11, nmax);

    for (size_t m = 0; m < 99; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdnpot
