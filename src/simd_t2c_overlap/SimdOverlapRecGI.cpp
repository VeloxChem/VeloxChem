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


#include "SimdOverlapRecGI.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"

#include "SimdOverlapVrrRecSD.hpp"
#include "SimdOverlapVrrRecSF.hpp"
#include "SimdOverlapVrrRecSG.hpp"
#include "SimdOverlapVrrRecSH.hpp"
#include "SimdOverlapVrrRecSI.hpp"
#include "SimdOverlapVrrRecSK.hpp"
#include "SimdOverlapVrrRecSL.hpp"
#include "SimdOverlapVrrRecSM.hpp"
#include "SimdOverlapVrrRecSN.hpp"
#include "SimdOverlapVrrRecSP.hpp"
#include "SimdOverlapVrrRecSS.hpp"
#include "SimdTransferDI.hpp"
#include "SimdTransferDK.hpp"
#include "SimdTransferDL.hpp"
#include "SimdTransferFI.hpp"
#include "SimdTransferFK.hpp"
#include "SimdTransferGI.hpp"
#include "SimdTransferPI.hpp"
#include "SimdTransferPK.hpp"
#include "SimdTransferPL.hpp"
#include "SimdTransferPM.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformI.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_gi_overlap(double               *values,
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
            false, std::string("compute_gi_overlap: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 2920, 289, 230, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 117 * nvalues, 0.0);

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

            const auto fb = a_exps[i] / p;

            simdfunc::compute_pb(buffer, coordinates, 0, ncols, fb);

            compute_prim_ss_overlap(buffer, coordinates, 3, ncols, fovl, mu);

            compute_prim_sp_overlap_0(buffer, 4, 0, 3, ncols);

            compute_prim_sd_overlap_0(buffer, 7, 0, 3, 4, ncols, p);

            compute_prim_sf_overlap_0(buffer, 13, 0, 4, 7, ncols, p);

            compute_prim_sg_overlap_0(buffer, 23, 0, 7, 13, ncols, p);

            compute_prim_sh_overlap_0(buffer, 38, 0, 13, 23, ncols, p);

            compute_prim_si_overlap_0(buffer, 59, 0, 23, 38, ncols, p);

            compute_prim_sk_overlap_0(buffer, 87, 0, 38, 59, ncols, p);

            compute_prim_sl_overlap_0(buffer, 123, 0, 59, 87, ncols, p);

            compute_prim_sm_overlap_0(buffer, 168, 0, 87, 123, ncols, p);

            compute_prim_sn_overlap_0(buffer, 223, 0, 123, 168, ncols, p);

            simdfunc::contract_primitives(buffer, 289, 59, 230, ncols);
        }
    }

    simdtrf::compute_hrr_pi(buffer, coordinates, 519, 289, 317, 1, nmax);

    simdtrf::compute_hrr_pk(buffer, coordinates, 603, 317, 353, 1, nmax);

    simdtrf::compute_hrr_pl(buffer, coordinates, 711, 353, 398, 1, nmax);

    simdtrf::compute_hrr_pm(buffer, coordinates, 846, 398, 453, 1, nmax);

    simdtrf::compute_hrr_di(buffer, coordinates, 1011, 519, 603, 1, nmax);

    simdtrf::compute_hrr_dk(buffer, coordinates, 1179, 603, 711, 1, nmax);

    simdtrf::compute_hrr_dl(buffer, coordinates, 1395, 711, 846, 1, nmax);

    simdtrf::compute_hrr_fi(buffer, coordinates, 1665, 1011, 1179, 1, nmax);

    simdtrf::compute_hrr_fk(buffer, coordinates, 1945, 1179, 1395, 1, nmax);

    simdtrf::compute_hrr_gi(buffer, coordinates, 2305, 1665, 1945, 1, nmax);

    simdtrf::transform_i_inner(buffer, 2725, 2305, 15, 1, nmax);

    simdtrf::transform_g_outer(values, nvalues, buffer, 2725, 13, nmax);

    for (size_t m = 0; m < 117; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdovl
