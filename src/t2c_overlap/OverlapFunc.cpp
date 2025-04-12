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

#include "OverlapFunc.hpp"

#include "OverlapRecDD.hpp"
#include "OverlapRecDF.hpp"
#include "OverlapRecDP.hpp"
#include "OverlapRecDS.hpp"
#include "OverlapRecFD.hpp"
#include "OverlapRecFF.hpp"
#include "OverlapRecFP.hpp"
#include "OverlapRecFS.hpp"
#include "OverlapRecPD.hpp"
#include "OverlapRecPF.hpp"
#include "OverlapRecPP.hpp"
#include "OverlapRecPS.hpp"
#include "OverlapRecSD.hpp"
#include "OverlapRecSF.hpp"
#include "OverlapRecSP.hpp"
#include "OverlapRecSS.hpp"

namespace ovlfunc {  // ovlfunc namespace

auto
compute(CSubMatrix* matrix, const CGtoBlock& gto_block, const int64_t angmom, const int64_t bra_first, const int64_t bra_last) -> void
{
    if (angmom == 0)
    {
        ovlrec::compOverlapSS(matrix, gto_block, bra_first, bra_last);

        return;
    }

    if (angmom == 1)
    {
        ovlrec::compOverlapPP(matrix, gto_block, bra_first, bra_last);

        return;
    }

    if (angmom == 2)
    {
        ovlrec::compOverlapDD(matrix, gto_block, bra_first, bra_last);

        return;
    }

    if (angmom == 3)
    {
        ovlrec::compOverlapFF(matrix, gto_block, bra_first, bra_last);

        return;
    }
}

auto
compute(CSubMatrix*      matrix,
        const CGtoBlock& bra_gto_block,
        const CGtoBlock& ket_gto_block,
        const int64_t    bra_angmom,
        const int64_t    ket_angmom,
        const bool       ang_order,
        const int64_t    bra_first,
        const int64_t    bra_last,
        const mat_t      mat_type) -> void
{
    if ((bra_angmom == 0) && (ket_angmom == 0))
    {
        ovlrec::compOverlapSS(matrix, bra_gto_block, ket_gto_block, bra_first, bra_last, mat_type);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 1))
    {
        ovlrec::compOverlapSP(matrix, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 0))
    {
        ovlrec::compOverlapPS(matrix, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 2))
    {
        ovlrec::compOverlapSD(matrix, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 0))
    {
        ovlrec::compOverlapDS(matrix, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 3))
    {
        ovlrec::compOverlapSF(matrix, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 0))
    {
        ovlrec::compOverlapFS(matrix, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 1))
    {
        ovlrec::compOverlapPP(matrix, bra_gto_block, ket_gto_block, bra_first, bra_last, mat_type);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 2))
    {
        ovlrec::compOverlapPD(matrix, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 1))
    {
        ovlrec::compOverlapDP(matrix, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 3))
    {
        ovlrec::compOverlapPF(matrix, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 1))
    {
        ovlrec::compOverlapFP(matrix, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 2))
    {
        ovlrec::compOverlapDD(matrix, bra_gto_block, ket_gto_block, bra_first, bra_last, mat_type);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 3))
    {
        ovlrec::compOverlapDF(matrix, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 2))
    {
        ovlrec::compOverlapFD(matrix, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 3))
    {
        ovlrec::compOverlapFF(matrix, bra_gto_block, ket_gto_block, bra_first, bra_last, mat_type);

        return;
    }
}

}  // namespace ovlfunc
