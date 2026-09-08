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



#include "SimdOverlapFunc.hpp"

#include <string>

#include "ErrorHandler.hpp"
#include "SimdOverlapRecSS.hpp"
#include "SimdOverlapRecSP.hpp"
#include "SimdOverlapRecSD.hpp"
#include "SimdOverlapRecPS.hpp"
#include "SimdOverlapRecPP.hpp"
#include "SimdOverlapRecPD.hpp"
#include "SimdOverlapRecDS.hpp"
#include "SimdOverlapRecDP.hpp"
#include "SimdOverlapRecDD.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_overlap(double               *values,
                const size_t          nvalues,
                const CBasisFunction &bra,
                const CBasisFunction &ket,
                const CSimdMatrix    &coordinates,
                const double          threshold) -> void
{
    const auto lbra = bra.get_angular_momentum();

    const auto lket = ket.get_angular_momentum();

    // NOTE: the two orders of a combination are separate kernels, as the recurrence
    // builds the angular momentum on ket side and transfers it to bra side, so the
    // work differs with the order even where the integrals do not.

    if ((lbra == 0) && (lket == 0))
    {
        compute_ss_overlap(values, nvalues, bra, ket, coordinates, threshold);

        return;
    }

    if ((lbra == 0) && (lket == 1))
    {
        compute_sp_overlap(values, nvalues, bra, ket, coordinates, threshold);

        return;
    }

    if ((lbra == 0) && (lket == 2))
    {
        compute_sd_overlap(values, nvalues, bra, ket, coordinates, threshold);

        return;
    }

    if ((lbra == 1) && (lket == 0))
    {
        compute_ps_overlap(values, nvalues, bra, ket, coordinates, threshold);

        return;
    }

    if ((lbra == 1) && (lket == 1))
    {
        compute_pp_overlap(values, nvalues, bra, ket, coordinates, threshold);

        return;
    }

    if ((lbra == 1) && (lket == 2))
    {
        compute_pd_overlap(values, nvalues, bra, ket, coordinates, threshold);

        return;
    }

    if ((lbra == 2) && (lket == 0))
    {
        compute_ds_overlap(values, nvalues, bra, ket, coordinates, threshold);

        return;
    }

    if ((lbra == 2) && (lket == 1))
    {
        compute_dp_overlap(values, nvalues, bra, ket, coordinates, threshold);

        return;
    }

    if ((lbra == 2) && (lket == 2))
    {
        compute_dd_overlap(values, nvalues, bra, ket, coordinates, threshold);

        return;
    }

    // NOTE: the kernels are generated up to angular momentum two, so a combination
    // above it stops rather than leaving the values of the sparsity pattern
    // unwritten, which is what a caller would otherwise read as integrals.

    errors::assertMsgCritical(false, std::string("SimdOverlapFunc.compute_overlap: Overlap integrals are not implemented"));
}

}  // namespace simdovl
