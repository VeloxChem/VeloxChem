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
#include "SimdOverlapRecSP.hpp"
#include "SimdOverlapRecPP.hpp"
#include "SimdOverlapRecPS.hpp"
#include "SimdOverlapRecSD.hpp"
#include "SimdOverlapRecDS.hpp"
#include "SimdOverlapRecSF.hpp"
#include "SimdOverlapRecFS.hpp"
#include "SimdOverlapRecSG.hpp"
#include "SimdOverlapRecGS.hpp"
#include "SimdOverlapRecSH.hpp"
#include "SimdOverlapRecHS.hpp"
#include "SimdOverlapRecSI.hpp"
#include "SimdOverlapRecIS.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_overlap(double                         *values,
                const size_t                    nvalues,
                const CBasisFunction           &bra,
                const CBasisFunction           &ket,
                const std::vector<CSimdMatrix> &harmonics,
                const CSimdMatrix              &coordinates,
                const double                    threshold) -> void
{
    const auto lbra = bra.get_angular_momentum();

    const auto lket = ket.get_angular_momentum();

    if ((lbra == 0) && (lket == 0))
    {
        compute_ss_overlap(values, nvalues, bra, ket, coordinates, threshold);

        return;
    }

    // NOTE: the solid harmonics of angular momentum l are the element of index
    // l - 1, as the harmonics of angular momentum zero are one for every atom
    // pair and are not stored.

    if (lbra == 0)
    {
        if (lket == 1)
        {
            compute_sp_overlap(values, nvalues, bra, ket, harmonics[0], coordinates, threshold);

            return;
        }

        if (lket == 2)
        {
            compute_sd_overlap(values, nvalues, bra, ket, harmonics[1], coordinates, threshold);

            return;
        }

        if (lket == 3)
        {
            compute_sf_overlap(values, nvalues, bra, ket, harmonics[2], coordinates, threshold);

            return;
        }

        if (lket == 4)
        {
            compute_sg_overlap(values, nvalues, bra, ket, harmonics[3], coordinates, threshold);

            return;
        }

        if (lket == 5)
        {
            compute_sh_overlap(values, nvalues, bra, ket, harmonics[4], coordinates, threshold);

            return;
        }

        if (lket == 6)
        {
            compute_si_overlap(values, nvalues, bra, ket, harmonics[5], coordinates, threshold);

            return;
        }

    }

    if ((lbra == 1) && (lket == 1))
    {
        compute_pp_overlap(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if (lket == 0)
    {
        if (lbra == 1)
        {
            compute_ps_overlap(values, nvalues, bra, ket, harmonics[0], coordinates, threshold);

            return;
        }

        if (lbra == 2)
        {
            compute_ds_overlap(values, nvalues, bra, ket, harmonics[1], coordinates, threshold);

            return;
        }

        if (lbra == 3)
        {
            compute_fs_overlap(values, nvalues, bra, ket, harmonics[2], coordinates, threshold);

            return;
        }

        if (lbra == 4)
        {
            compute_gs_overlap(values, nvalues, bra, ket, harmonics[3], coordinates, threshold);

            return;
        }

        if (lbra == 5)
        {
            compute_hs_overlap(values, nvalues, bra, ket, harmonics[4], coordinates, threshold);

            return;
        }

        if (lbra == 6)
        {
            compute_is_overlap(values, nvalues, bra, ket, harmonics[5], coordinates, threshold);

            return;
        }

    }

    errors::assertMsgCritical(
        false, std::string("SimdOverlapFunc.compute_overlap: Overlap integrals of the requested angular momenta are not implemented"));
}

}  // namespace simdovl
