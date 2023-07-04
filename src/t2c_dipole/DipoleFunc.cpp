#include "DipoleFunc.hpp"

#include "DipoleRecDD.hpp"
#include "DipoleRecDF.hpp"
#include "DipoleRecDG.hpp"
#include "DipoleRecDP.hpp"
#include "DipoleRecDS.hpp"
#include "DipoleRecFD.hpp"
#include "DipoleRecFF.hpp"
#include "DipoleRecFG.hpp"
#include "DipoleRecFP.hpp"
#include "DipoleRecFS.hpp"
#include "DipoleRecGD.hpp"
#include "DipoleRecGF.hpp"
#include "DipoleRecGG.hpp"
#include "DipoleRecGP.hpp"
#include "DipoleRecGS.hpp"
#include "DipoleRecPD.hpp"
#include "DipoleRecPF.hpp"
#include "DipoleRecPG.hpp"
#include "DipoleRecPP.hpp"
#include "DipoleRecPS.hpp"
#include "DipoleRecSD.hpp"
#include "DipoleRecSF.hpp"
#include "DipoleRecSG.hpp"
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
        mpol::compDipoleSS(matrix_x, matrix_y, matrix_z, point, gto_block, bra_first, bra_last);

        return;
    }

    if (angmom == 1)
    {
        mpol::compDipolePP(matrix_x, matrix_y, matrix_z, point, gto_block, bra_first, bra_last);

        return;
    }

    if (angmom == 2)
    {
        mpol::compDipoleDD(matrix_x, matrix_y, matrix_z, point, gto_block, bra_first, bra_last);

        return;
    }

    if (angmom == 3)
    {
        mpol::compDipoleFF(matrix_x, matrix_y, matrix_z, point, gto_block, bra_first, bra_last);

        return;
    }

    if (angmom == 4)
    {
        mpol::compDipoleGG(matrix_x, matrix_y, matrix_z, point, gto_block, bra_first, bra_last);

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
        mpol::compDipoleSS(matrix_x, matrix_y, matrix_z, point, bra_gto_block, ket_gto_block, bra_first, bra_last, mat_type);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 1))
    {
        mpol::compDipoleSP(matrix_x, matrix_y, matrix_z, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 0))
    {
        mpol::compDipolePS(matrix_x, matrix_y, matrix_z, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 2))
    {
        mpol::compDipoleSD(matrix_x, matrix_y, matrix_z, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 0))
    {
        mpol::compDipoleDS(matrix_x, matrix_y, matrix_z, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 3))
    {
        mpol::compDipoleSF(matrix_x, matrix_y, matrix_z, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 0))
    {
        mpol::compDipoleFS(matrix_x, matrix_y, matrix_z, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 4))
    {
        mpol::compDipoleSG(matrix_x, matrix_y, matrix_z, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 0))
    {
        mpol::compDipoleGS(matrix_x, matrix_y, matrix_z, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 1))
    {
        mpol::compDipolePP(matrix_x, matrix_y, matrix_z, point, bra_gto_block, ket_gto_block, bra_first, bra_last, mat_type);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 2))
    {
        mpol::compDipolePD(matrix_x, matrix_y, matrix_z, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 1))
    {
        mpol::compDipoleDP(matrix_x, matrix_y, matrix_z, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 3))
    {
        mpol::compDipolePF(matrix_x, matrix_y, matrix_z, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 1))
    {
        mpol::compDipoleFP(matrix_x, matrix_y, matrix_z, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 4))
    {
        mpol::compDipolePG(matrix_x, matrix_y, matrix_z, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 1))
    {
        mpol::compDipoleGP(matrix_x, matrix_y, matrix_z, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 2))
    {
        mpol::compDipoleDD(matrix_x, matrix_y, matrix_z, point, bra_gto_block, ket_gto_block, bra_first, bra_last, mat_type);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 3))
    {
        mpol::compDipoleDF(matrix_x, matrix_y, matrix_z, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 2))
    {
        mpol::compDipoleFD(matrix_x, matrix_y, matrix_z, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 4))
    {
        mpol::compDipoleDG(matrix_x, matrix_y, matrix_z, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 2))
    {
        mpol::compDipoleGD(matrix_x, matrix_y, matrix_z, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 3))
    {
        mpol::compDipoleFF(matrix_x, matrix_y, matrix_z, point, bra_gto_block, ket_gto_block, bra_first, bra_last, mat_type);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 4))
    {
        mpol::compDipoleFG(matrix_x, matrix_y, matrix_z, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 3))
    {
        mpol::compDipoleGF(matrix_x, matrix_y, matrix_z, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 4))
    {
        mpol::compDipoleGG(matrix_x, matrix_y, matrix_z, point, bra_gto_block, ket_gto_block, bra_first, bra_last, mat_type);

        return;
    }
}

}  // namespace dipfunc
