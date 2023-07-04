#include "DipoleFunc.hpp"

#include "DipoleRecDD.hpp"
#include "DipoleRecDF.hpp"
#include "DipoleRecDP.hpp"
#include "DipoleRecDS.hpp"
#include "DipoleRecFD.hpp"
#include "DipoleRecFF.hpp"
#include "DipoleRecFP.hpp"
#include "DipoleRecFS.hpp"
#include "DipoleRecPD.hpp"
#include "DipoleRecPF.hpp"
#include "DipoleRecPP.hpp"
#include "DipoleRecPS.hpp"
#include "DipoleRecSD.hpp"
#include "DipoleRecSF.hpp"
#include "DipoleRecSP.hpp"
#include "DipoleRecSS.hpp"

namespace dipfunc {  // dipfunc namespace

auto
compute(CSubMatrix*      matrix_x,
        CSubMatrix*      matrix_y,
        CSubMatrix*      matrix_z,
        const TPoint3D&  point,
        const CGtoBlock& gto_block,
        const int64_t    angmom,
        const int64_t    bra_first,
        const int64_t    bra_last) -> void
{
    if (angmom == 0)
    {
        diprec::compDipoleSS(matrix_x, matrix_y, matrix_z, point, gto_block, bra_first, bra_last);

        return;
    }

    if (angmom == 1)
    {
        diprec::compDipolePP(matrix_x, matrix_y, matrix_z, point, gto_block, bra_first, bra_last);

        return;
    }

    if (angmom == 2)
    {
        diprec::compDipoleDD(matrix_x, matrix_y, matrix_z, point, gto_block, bra_first, bra_last);

        return;
    }

    if (angmom == 3)
    {
        diprec::compDipoleFF(matrix_x, matrix_y, matrix_z, point, gto_block, bra_first, bra_last);

        return;
    }
}

auto
compute(CSubMatrix*      matrix_x,
        CSubMatrix*      matrix_y,
        CSubMatrix*      matrix_z,
        const TPoint3D&  point,
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
        diprec::compDipoleSS(matrix_x, matrix_y, matrix_z, point, bra_gto_block, ket_gto_block, bra_first, bra_last, mat_type);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 1))
    {
        diprec::compDipoleSP(matrix_x, matrix_y, matrix_z, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 0))
    {
        diprec::compDipolePS(matrix_x, matrix_y, matrix_z, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 2))
    {
        diprec::compDipoleSD(matrix_x, matrix_y, matrix_z, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 0))
    {
        diprec::compDipoleDS(matrix_x, matrix_y, matrix_z, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 3))
    {
        diprec::compDipoleSF(matrix_x, matrix_y, matrix_z, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 0))
    {
        diprec::compDipoleFS(matrix_x, matrix_y, matrix_z, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 1))
    {
        diprec::compDipolePP(matrix_x, matrix_y, matrix_z, point, bra_gto_block, ket_gto_block, bra_first, bra_last, mat_type);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 2))
    {
        diprec::compDipolePD(matrix_x, matrix_y, matrix_z, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 1))
    {
        diprec::compDipoleDP(matrix_x, matrix_y, matrix_z, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 3))
    {
        diprec::compDipolePF(matrix_x, matrix_y, matrix_z, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 1))
    {
        diprec::compDipoleFP(matrix_x, matrix_y, matrix_z, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 2))
    {
        diprec::compDipoleDD(matrix_x, matrix_y, matrix_z, point, bra_gto_block, ket_gto_block, bra_first, bra_last, mat_type);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 3))
    {
        diprec::compDipoleDF(matrix_x, matrix_y, matrix_z, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 2))
    {
        diprec::compDipoleFD(matrix_x, matrix_y, matrix_z, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 3))
    {
        diprec::compDipoleFF(matrix_x, matrix_y, matrix_z, point, bra_gto_block, ket_gto_block, bra_first, bra_last, mat_type);

        return;
    }
}

}  // namespace dipfunc
