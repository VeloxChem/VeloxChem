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


#include "SimdNuclearPotentialRecIP.hpp"

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

#include "SimdNuclearPotentialVrrRecDS.hpp"
#include "SimdNuclearPotentialVrrRecFS.hpp"
#include "SimdNuclearPotentialVrrRecGS.hpp"
#include "SimdNuclearPotentialVrrRecHS.hpp"
#include "SimdNuclearPotentialVrrRecIS.hpp"
#include "SimdNuclearPotentialVrrRecKS.hpp"
#include "SimdNuclearPotentialVrrRecPS.hpp"
#include "SimdTransferIP.hpp"
#include "SimdTransformI.hpp"
#include "SimdTransformP.hpp"

namespace simdnpot {  // simdnpot namespace

auto
compute_ip_nuclear_potential(double                    *values,
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
            false, std::string("compute_ip_nuclear_potential: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    errors::assertMsgCritical(
        points.size() == 3 * charges.size(),
        std::string("compute_ip_nuclear_potential: Expecting three coordinates for each charge"));

    if (charges.empty())
    {
        std::fill(values, values + 39 * nvalues, 0.0);

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

    const auto nmax = simdfunc::prepare_buffer(buffer, 569, 337, 64, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 39 * nvalues, 0.0);

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

            const auto fa = -b_exps[j] / p;

            const auto fc = b_exps[j] / p;

            simdfunc::compute_pa(buffer, coordinates, 0, ncols, fa);

            for (size_t ic = 0; ic < charges.size(); ic++)
            {
                const auto fz = fnpot * charges[ic];

                simdfunc::compute_pc(buffer, coordinates, 3, points, ic, ncols, fc);

                simdfunc::compute_full_npot_boys_function(buffer, coordinates, 6, 3, 7, ncols,
                                                          fz, mu, p);

                compute_prim_ps_nuclear_potential_0(buffer, 15, 0, 3, 7, 8, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 18, 0, 3, 8, 9, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 21, 0, 3, 9, 10, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 24, 0, 3, 10, 11, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 27, 0, 3, 11, 12, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 30, 0, 3, 12, 13, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 33, 0, 3, 13, 14, ncols);

                compute_prim_ds_nuclear_potential_0(buffer, 36, 0, 3, 7, 8, 15, 18, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 42, 0, 3, 8, 9, 18, 21, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 48, 0, 3, 9, 10, 21, 24, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 54, 0, 3, 10, 11, 24, 27, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 60, 0, 3, 11, 12, 27, 30, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 66, 0, 3, 12, 13, 30, 33, ncols, p);

                compute_prim_fs_nuclear_potential_0(buffer, 72, 0, 3, 15, 18, 36, 42, ncols, p);

                compute_prim_fs_nuclear_potential_0(buffer, 82, 0, 3, 18, 21, 42, 48, ncols, p);

                compute_prim_fs_nuclear_potential_0(buffer, 92, 0, 3, 21, 24, 48, 54, ncols, p);

                compute_prim_fs_nuclear_potential_0(buffer, 102, 0, 3, 24, 27, 54, 60, ncols,
                                                    p);

                compute_prim_fs_nuclear_potential_0(buffer, 112, 0, 3, 27, 30, 60, 66, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 122, 0, 3, 36, 42, 72, 82, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 137, 0, 3, 42, 48, 82, 92, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 152, 0, 3, 48, 54, 92, 102, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 167, 0, 3, 54, 60, 102, 112, ncols,
                                                    p);

                compute_prim_hs_nuclear_potential_0(buffer, 182, 0, 3, 72, 82, 122, 137, ncols,
                                                    p);

                compute_prim_hs_nuclear_potential_0(buffer, 203, 0, 3, 82, 92, 137, 152, ncols,
                                                    p);

                compute_prim_hs_nuclear_potential_0(buffer, 224, 0, 3, 92, 102, 152, 167, ncols,
                                                    p);

                compute_prim_is_nuclear_potential_0(buffer, 245, 0, 3, 122, 137, 182, 203, ncols,
                                                    p);

                compute_prim_is_nuclear_potential_0(buffer, 273, 0, 3, 137, 152, 203, 224, ncols,
                                                    p);

                compute_prim_ks_nuclear_potential_0(buffer, 301, 0, 3, 182, 203, 245, 273, ncols,
                                                    p);

                simdfunc::contract_primitives(buffer, 337, 245, 28, ncols);

                simdfunc::contract_primitives(buffer, 365, 301, 36, ncols);
            }
        }
    }

    simdtrf::compute_hrr_ip(buffer, coordinates, 401, 337, 365, 1, nmax);

    simdtrf::transform_p_inner(buffer, 485, 401, 28, 1, nmax);

    simdtrf::transform_i_outer(values, nvalues, buffer, 485, 3, nmax);

    for (size_t m = 0; m < 39; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdnpot
