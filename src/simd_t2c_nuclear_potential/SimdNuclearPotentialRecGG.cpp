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

    const auto nmax = simdfunc::prepare_buffer(buffer, 2051, 502, 145, dimensions);

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

            for (size_t ic = 0; ic < charges.size(); ic++)
            {
                const auto fz = fnpot * charges[ic];

                simdfunc::compute_pc(buffer, coordinates, 3, points, ic, ncols, fc);

                simdfunc::compute_full_npot_boys_function(buffer, coordinates, 6, 3, 8, ncols,
                                                          fz, mu, p);

                compute_prim_sp_nuclear_potential_0(buffer, 16, 0, 3, 7, 8, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 19, 0, 3, 8, 9, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 22, 0, 3, 9, 10, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 25, 0, 3, 10, 11, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 28, 0, 3, 11, 12, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 31, 0, 3, 12, 13, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 34, 0, 3, 13, 14, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 37, 0, 3, 14, 15, ncols);

                compute_prim_sd_nuclear_potential_0(buffer, 40, 0, 3, 7, 16, 8, 19, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 46, 0, 3, 8, 19, 9, 22, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 52, 0, 3, 9, 22, 10, 25, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 58, 0, 3, 10, 25, 11, 28, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 64, 0, 3, 11, 28, 12, 31, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 70, 0, 3, 12, 31, 13, 34, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 76, 0, 3, 13, 34, 14, 37, ncols, p);

                compute_prim_sf_nuclear_potential_0(buffer, 82, 0, 3, 16, 40, 19, 46, ncols, p);

                compute_prim_sf_nuclear_potential_0(buffer, 92, 0, 3, 19, 46, 22, 52, ncols, p);

                compute_prim_sf_nuclear_potential_0(buffer, 102, 0, 3, 22, 52, 25, 58, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 112, 0, 3, 25, 58, 28, 64, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 122, 0, 3, 28, 64, 31, 70, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 132, 0, 3, 31, 70, 34, 76, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 142, 0, 3, 40, 82, 46, 92, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 157, 0, 3, 46, 92, 52, 102, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 172, 0, 3, 52, 102, 58, 112, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 187, 0, 3, 58, 112, 64, 122, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 202, 0, 3, 64, 122, 70, 132, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 217, 0, 3, 82, 142, 92, 157, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 238, 0, 3, 92, 157, 102, 172, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 259, 0, 3, 102, 172, 112, 187, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 280, 0, 3, 112, 187, 122, 202, ncols,
                                                    p);

                compute_prim_si_nuclear_potential_0(buffer, 301, 0, 3, 142, 217, 157, 238, ncols,
                                                    p);

                compute_prim_si_nuclear_potential_0(buffer, 329, 0, 3, 157, 238, 172, 259, ncols,
                                                    p);

                compute_prim_si_nuclear_potential_0(buffer, 357, 0, 3, 172, 259, 187, 280, ncols,
                                                    p);

                compute_prim_sk_nuclear_potential_0(buffer, 385, 0, 3, 217, 301, 238, 329, ncols,
                                                    p);

                compute_prim_sk_nuclear_potential_0(buffer, 421, 0, 3, 238, 329, 259, 357, ncols,
                                                    p);

                compute_prim_sl_nuclear_potential_0(buffer, 457, 0, 3, 301, 385, 329, 421, ncols,
                                                    p);

                simdfunc::contract_primitives(buffer, 502, 142, 15, ncols);

                simdfunc::contract_primitives(buffer, 517, 217, 21, ncols);

                simdfunc::contract_primitives(buffer, 538, 301, 28, ncols);

                simdfunc::contract_primitives(buffer, 566, 385, 36, ncols);

                simdfunc::contract_primitives(buffer, 602, 457, 45, ncols);
            }
        }
    }

    simdtrf::compute_hrr_pg(buffer, coordinates, 647, 502, 517, 1, nmax);

    simdtrf::compute_hrr_ph(buffer, coordinates, 692, 517, 538, 1, nmax);

    simdtrf::compute_hrr_pi(buffer, coordinates, 755, 538, 566, 1, nmax);

    simdtrf::compute_hrr_pk(buffer, coordinates, 839, 566, 602, 1, nmax);

    simdtrf::compute_hrr_dg(buffer, coordinates, 947, 647, 692, 1, nmax);

    simdtrf::compute_hrr_dh(buffer, coordinates, 1037, 692, 755, 1, nmax);

    simdtrf::compute_hrr_di(buffer, coordinates, 1163, 755, 839, 1, nmax);

    simdtrf::compute_hrr_fg(buffer, coordinates, 1331, 947, 1037, 1, nmax);

    simdtrf::compute_hrr_fh(buffer, coordinates, 1481, 1037, 1163, 1, nmax);

    simdtrf::compute_hrr_gg(buffer, coordinates, 1691, 1331, 1481, 1, nmax);

    simdtrf::transform_g_inner(buffer, 1916, 1691, 15, 1, nmax);

    simdtrf::transform_g_outer(values, nvalues, buffer, 1916, 9, nmax);

    for (size_t m = 0; m < 81; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdnpot
