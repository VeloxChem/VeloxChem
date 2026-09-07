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

#include <algorithm>
#include <string>

#include "ErrorHandler.hpp"
#include "SimdOverlapRecSS.hpp"
#include "SimdOverlapRecSLP.hpp"
#include "SimdOverlapRecSLD.hpp"
#include "SimdOverlapRecSLF.hpp"
#include "SimdOverlapRecSLG.hpp"
#include "SimdOverlapRecSLH.hpp"
#include "SimdOverlapRecSLI.hpp"

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

    if ((lbra == 0) && (lket == 0))
    {
        compute_ss_overlap(values, nvalues, bra, ket, coordinates, threshold);

        return;
    }

    // NOTE: a combination of one basis function of zero angular momentum and one of
    // higher angular momentum is computed by a single kernel in either order, as the
    // harmonic is the same polynomial of the vector between the atoms either way and
    // the orders differ only in the prefactor.

    const auto lmin = std::min(lbra, lket);

    const auto lmax = std::max(lbra, lket);

    if (lmin == 0)
    {
        if (lmax == 1)
        {
            compute_slp_overlap(values, nvalues, bra, ket, coordinates, threshold);

            return;
        }

        if (lmax == 2)
        {
            compute_sld_overlap(values, nvalues, bra, ket, coordinates, threshold);

            return;
        }

        if (lmax == 3)
        {
            compute_slf_overlap(values, nvalues, bra, ket, coordinates, threshold);

            return;
        }

        if (lmax == 4)
        {
            compute_slg_overlap(values, nvalues, bra, ket, coordinates, threshold);

            return;
        }

        if (lmax == 5)
        {
            compute_slh_overlap(values, nvalues, bra, ket, coordinates, threshold);

            return;
        }

        if (lmax == 6)
        {
            compute_sli_overlap(values, nvalues, bra, ket, coordinates, threshold);

            return;
        }

    }

    // NOTE: the kernels of the remaining combinations of basis functions are being
    // rewritten, so they are not computed here. A combination stops rather than
    // leaving the values of the sparsity pattern unwritten, which is what a caller
    // would otherwise read as integrals.

    errors::assertMsgCritical(false, std::string("SimdOverlapFunc.compute_overlap: Overlap integrals are not implemented"));
}

}  // namespace simdovl
