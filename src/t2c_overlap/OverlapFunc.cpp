//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

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
