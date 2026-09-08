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


#include "SimdOverlapRecII.hpp"

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
#include "SimdOverlapVrrRecSO.hpp"
#include "SimdOverlapVrrRecSP.hpp"
#include "SimdOverlapVrrRecSQ.hpp"
#include "SimdOverlapVrrRecSS.hpp"
#include "SimdTransferDI.hpp"
#include "SimdTransferDK.hpp"
#include "SimdTransferDL.hpp"
#include "SimdTransferDM.hpp"
#include "SimdTransferDN.hpp"
#include "SimdTransferFI.hpp"
#include "SimdTransferFK.hpp"
#include "SimdTransferFL.hpp"
#include "SimdTransferFM.hpp"
#include "SimdTransferGI.hpp"
#include "SimdTransferGK.hpp"
#include "SimdTransferGL.hpp"
#include "SimdTransferHI.hpp"
#include "SimdTransferHK.hpp"
#include "SimdTransferII.hpp"
#include "SimdTransferPI.hpp"
#include "SimdTransferPK.hpp"
#include "SimdTransferPL.hpp"
#include "SimdTransferPM.hpp"
#include "SimdTransferPN.hpp"
#include "SimdTransferPO.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_ii_overlap(double               *values,
                        const size_t          nvalues,
                        const CBasisFunction &bra,
                        const CBasisFunction &ket,
                        const CSimdMatrix    &coordinates,
                        const double          threshold) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_ii_overlap: Number of values exceeds number of atom pairs"));
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

    auto buffer = simdfunc::make_primitive_buffer(dimensions, 7761);

    if (buffer.number_of_columns() == 0)
    {
        std::fill(values, values + 169 * nvalues, 0.0);

        return;
    }

    const auto nmax = buffer.number_of_columns();

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

            compute_prim_sd_overlap_1(buffer, 7, 0, 3, 4, ncols, p);

            compute_prim_sf_overlap_1(buffer, 10, 0, 4, 7, ncols, p);

            compute_prim_sg_overlap_1(buffer, 16, 0, 7, 10, ncols, p);

            compute_prim_sh_overlap_0(buffer, 25, 0, 10, 16, ncols, p);

            compute_prim_si_overlap_1(buffer, 40, 0, 16, 25, ncols, p);

            compute_prim_sk_overlap_0(buffer, 68, 0, 25, 40, ncols, p);

            compute_prim_sl_overlap_0(buffer, 104, 0, 40, 68, ncols, p);

            compute_prim_sm_overlap_0(buffer, 149, 0, 68, 104, ncols, p);

            compute_prim_sn_overlap_0(buffer, 204, 0, 104, 149, ncols, p);

            compute_prim_so_overlap_0(buffer, 270, 0, 149, 204, ncols, p);

            compute_prim_sq_overlap_0(buffer, 348, 0, 204, 270, ncols, p);

            simdfunc::contract_primitives(buffer, 439, 40, 399, ncols);
        }
    }

    compute_hrr_pi(buffer, coordinates, 838, 439, 467, nmax);

    compute_hrr_pk(buffer, coordinates, 922, 467, 503, nmax);

    compute_hrr_pl(buffer, coordinates, 1030, 503, 548, nmax);

    compute_hrr_pm(buffer, coordinates, 1165, 548, 603, nmax);

    compute_hrr_pn(buffer, coordinates, 1330, 603, 669, nmax);

    compute_hrr_po(buffer, coordinates, 1528, 669, 747, nmax);

    compute_hrr_di(buffer, coordinates, 1762, 838, 922, nmax);

    compute_hrr_dk(buffer, coordinates, 1930, 922, 1030, nmax);

    compute_hrr_dl(buffer, coordinates, 2146, 1030, 1165, nmax);

    compute_hrr_dm(buffer, coordinates, 2416, 1165, 1330, nmax);

    compute_hrr_dn(buffer, coordinates, 2746, 1330, 1528, nmax);

    compute_hrr_fi(buffer, coordinates, 3142, 1762, 1930, nmax);

    compute_hrr_fk(buffer, coordinates, 3422, 1930, 2146, nmax);

    compute_hrr_fl(buffer, coordinates, 3782, 2146, 2416, nmax);

    compute_hrr_fm(buffer, coordinates, 4232, 2416, 2746, nmax);

    compute_hrr_gi(buffer, coordinates, 4782, 3142, 3422, nmax);

    compute_hrr_gk(buffer, coordinates, 5202, 3422, 3782, nmax);

    compute_hrr_gl(buffer, coordinates, 5742, 3782, 4232, nmax);

    compute_hrr_hi(buffer, coordinates, 6417, 4782, 5202, nmax);

    compute_hrr_hk(buffer, coordinates, 7005, 5202, 5742, nmax);

    compute_hrr_ii_sph_tri(values, nvalues, buffer, coordinates, 6417, 7005, nmax);

    for (size_t m = 0; m < 169; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdovl
