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


#include "SimdNuclearPotentialRecID.hpp"

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
#include "SimdNuclearPotentialVrrRecLS.hpp"
#include "SimdNuclearPotentialVrrRecPS.hpp"
#include "SimdTransferID.hpp"
#include "SimdTransferIP.hpp"
#include "SimdTransferKP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformI.hpp"

namespace simdnpot {  // simdnpot namespace

auto
compute_id_nuclear_potential(double                    *values,
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
            false, std::string("compute_id_nuclear_potential: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    errors::assertMsgCritical(
        points.size() == 3 * charges.size(),
        std::string("compute_id_nuclear_potential: Expecting three coordinates for each charge"));

    if (charges.empty())
    {
        std::fill(values, values + 65 * nvalues, 0.0);

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

    const auto nmax = simdfunc::prepare_buffer(buffer, 1112, 503, 109, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 65 * nvalues, 0.0);

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

            simdfunc::compute_pair_exponent(buffer, coordinates, 6, ncols, mu);

            for (size_t ic = 0; ic < charges.size(); ic++)
            {
                const auto fz = fnpot * charges[ic];

                simdfunc::compute_pc(buffer, coordinates, 3, points, ic, ncols, fc);

                simdfunc::compute_full_npot_boys_function(buffer, coordinates, 7, 3, 8, ncols,
                                                          fz, 6, p);

                compute_prim_ps_nuclear_potential_0(buffer, 17, 0, 3, 8, 9, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 20, 0, 3, 9, 10, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 23, 0, 3, 10, 11, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 26, 0, 3, 11, 12, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 29, 0, 3, 12, 13, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 32, 0, 3, 13, 14, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 35, 0, 3, 14, 15, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 38, 0, 3, 15, 16, ncols);

                compute_prim_ds_nuclear_potential_0(buffer, 41, 0, 3, 8, 9, 17, 20, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 47, 0, 3, 9, 10, 20, 23, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 53, 0, 3, 10, 11, 23, 26, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 59, 0, 3, 11, 12, 26, 29, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 65, 0, 3, 12, 13, 29, 32, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 71, 0, 3, 13, 14, 32, 35, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 77, 0, 3, 14, 15, 35, 38, ncols, p);

                compute_prim_fs_nuclear_potential_0(buffer, 83, 0, 3, 17, 20, 41, 47, ncols, p);

                compute_prim_fs_nuclear_potential_0(buffer, 93, 0, 3, 20, 23, 47, 53, ncols, p);

                compute_prim_fs_nuclear_potential_0(buffer, 103, 0, 3, 23, 26, 53, 59, ncols,
                                                    p);

                compute_prim_fs_nuclear_potential_0(buffer, 113, 0, 3, 26, 29, 59, 65, ncols,
                                                    p);

                compute_prim_fs_nuclear_potential_0(buffer, 123, 0, 3, 29, 32, 65, 71, ncols,
                                                    p);

                compute_prim_fs_nuclear_potential_0(buffer, 133, 0, 3, 32, 35, 71, 77, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 143, 0, 3, 41, 47, 83, 93, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 158, 0, 3, 47, 53, 93, 103, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 173, 0, 3, 53, 59, 103, 113, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 188, 0, 3, 59, 65, 113, 123, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 203, 0, 3, 65, 71, 123, 133, ncols,
                                                    p);

                compute_prim_hs_nuclear_potential_0(buffer, 218, 0, 3, 83, 93, 143, 158, ncols,
                                                    p);

                compute_prim_hs_nuclear_potential_0(buffer, 239, 0, 3, 93, 103, 158, 173, ncols,
                                                    p);

                compute_prim_hs_nuclear_potential_0(buffer, 260, 0, 3, 103, 113, 173, 188, ncols,
                                                    p);

                compute_prim_hs_nuclear_potential_0(buffer, 281, 0, 3, 113, 123, 188, 203, ncols,
                                                    p);

                compute_prim_is_nuclear_potential_0(buffer, 302, 0, 3, 143, 158, 218, 239, ncols,
                                                    p);

                compute_prim_is_nuclear_potential_0(buffer, 330, 0, 3, 158, 173, 239, 260, ncols,
                                                    p);

                compute_prim_is_nuclear_potential_0(buffer, 358, 0, 3, 173, 188, 260, 281, ncols,
                                                    p);

                compute_prim_ks_nuclear_potential_0(buffer, 386, 0, 3, 218, 239, 302, 330, ncols,
                                                    p);

                compute_prim_ks_nuclear_potential_0(buffer, 422, 0, 3, 239, 260, 330, 358, ncols,
                                                    p);

                compute_prim_ls_nuclear_potential_0(buffer, 458, 0, 3, 302, 330, 386, 422, ncols,
                                                    p);

                simdfunc::contract_primitives(buffer, 503, 302, 28, ncols);

                simdfunc::contract_primitives(buffer, 531, 386, 36, ncols);

                simdfunc::contract_primitives(buffer, 567, 458, 45, ncols);
            }
        }
    }

    simdtrf::compute_hrr_ip(buffer, coordinates, 612, 503, 531, 1, nmax);

    simdtrf::compute_hrr_kp(buffer, coordinates, 696, 531, 567, 1, nmax);

    simdtrf::compute_hrr_id(buffer, coordinates, 804, 612, 696, 1, nmax);

    simdtrf::transform_d_inner(buffer, 972, 804, 28, 1, nmax);

    simdtrf::transform_i_outer(values, nvalues, buffer, 972, 5, nmax);

    for (size_t m = 0; m < 65; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdnpot
