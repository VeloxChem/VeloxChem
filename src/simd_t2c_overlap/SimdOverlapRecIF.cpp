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


#include "SimdOverlapRecIF.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"

#include "SimdOverlapVrrRecDS.hpp"
#include "SimdOverlapVrrRecFS.hpp"
#include "SimdOverlapVrrRecGS.hpp"
#include "SimdOverlapVrrRecHS.hpp"
#include "SimdOverlapVrrRecIS.hpp"
#include "SimdOverlapVrrRecKS.hpp"
#include "SimdOverlapVrrRecLS.hpp"
#include "SimdOverlapVrrRecMS.hpp"
#include "SimdOverlapVrrRecPS.hpp"
#include "SimdOverlapVrrRecSS.hpp"
#include "SimdTransferID.hpp"
#include "SimdTransferIF.hpp"
#include "SimdTransferIP.hpp"
#include "SimdTransferKD.hpp"
#include "SimdTransferKP.hpp"
#include "SimdTransferLP.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformI.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_if_overlap(double               *values,
                   const size_t          nvalues,
                   const CBasisFunction &bra,
                   const CBasisFunction &ket,
                   const CSimdMatrix    &coordinates,
                   CSimdMatrix          &buffer,
                   const double          threshold) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_if_overlap: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nprims = nprim_a * nprim_b;

    // NOTE: the pairs of primitives are screened with the threshold of the
    // integrals divided by their number, as every integral is a sum over them
    // and the error of a sum is bounded by the number of its terms.

    const auto dimensions = simdfunc::make_column_dimensions(
        bra, ket, nvalues, coordinates, screenfunc::two_center_overlap_primitive_bound, threshold / static_cast<double>(nprims));

    const auto nmax = simdfunc::prepare_buffer(buffer, 1574, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 91 * nvalues, 0.0);

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

            const auto fpi = mathconst::pi_value() / p;

            const auto fovl = a_norms[i] * b_norms[j] * fpi * std::sqrt(fpi);

            const auto fa = -b_exps[j] / p;

            simdfunc::compute_pa(buffer, coordinates, 0, ncols, fa);

            compute_prim_ss_overlap(buffer, coordinates, 3, ncols, fovl, mu);

            compute_prim_ps_overlap_0(buffer, 4, 0, 3, ncols);

            compute_prim_ds_overlap_0(buffer, 7, 0, 3, 4, ncols, p);

            compute_prim_fs_overlap_0(buffer, 13, 0, 4, 7, ncols, p);

            compute_prim_gs_overlap_0(buffer, 23, 0, 7, 13, ncols, p);

            compute_prim_hs_overlap_0(buffer, 38, 0, 13, 23, ncols, p);

            compute_prim_is_overlap_0(buffer, 59, 0, 23, 38, ncols, p);

            compute_prim_ks_overlap_0(buffer, 87, 0, 38, 59, ncols, p);

            compute_prim_ls_overlap_0(buffer, 123, 0, 59, 87, ncols, p);

            compute_prim_ms_overlap_0(buffer, 168, 0, 87, 123, ncols, p);

            simdfunc::contract_primitives(buffer, 223, 59, 164, ncols);
        }
    }

    simdtrf::compute_hrr_ip(buffer, coordinates, 387, 223, 251, nmax);

    simdtrf::compute_hrr_kp(buffer, coordinates, 471, 251, 287, nmax);

    simdtrf::compute_hrr_lp(buffer, coordinates, 579, 287, 332, nmax);

    simdtrf::compute_hrr_id(buffer, coordinates, 714, 387, 471, nmax);

    simdtrf::compute_hrr_kd(buffer, coordinates, 882, 471, 579, nmax);

    simdtrf::compute_hrr_if(buffer, coordinates, 1098, 714, 882, nmax);

    simdtrf::transform_f_inner(buffer, 1378, 1098, 28, nmax);

    simdtrf::transform_i_outer(values, nvalues, buffer, 1378, 7, nmax);

    for (size_t m = 0; m < 91; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdovl
