#include "QuadrupoleFunc.hpp"

#include "QuadrupoleRecDD.hpp"
#include "QuadrupoleRecDF.hpp"
#include "QuadrupoleRecDP.hpp"
#include "QuadrupoleRecDS.hpp"
#include "QuadrupoleRecFD.hpp"
#include "QuadrupoleRecFF.hpp"
#include "QuadrupoleRecFP.hpp"
#include "QuadrupoleRecFS.hpp"
#include "QuadrupoleRecPD.hpp"
#include "QuadrupoleRecPF.hpp"
#include "QuadrupoleRecPP.hpp"
#include "QuadrupoleRecPS.hpp"
#include "QuadrupoleRecSD.hpp"
#include "QuadrupoleRecSF.hpp"
#include "QuadrupoleRecSP.hpp"
#include "QuadrupoleRecSS.hpp"

namespace quadfunc {  // quadfunc namespace

auto
compute(CSubMatrix*      matrix_xx,
        CSubMatrix*      matrix_xy,
        CSubMatrix*      matrix_xz,
        CSubMatrix*      matrix_yy,
        CSubMatrix*      matrix_yz,
        CSubMatrix*      matrix_zz,
        const TPoint3D&  point,
        const CGtoBlock& gto_block,
        const int64_t    angmom,
        const int64_t    bra_first,
        const int64_t    bra_last) -> void
{
    if (angmom == 0)
    {
        quadrec::compQuadrupoleSS(matrix_xx, matrix_xy, matrix_xz, matrix_yy, matrix_yz, matrix_zz, point, gto_block, bra_first, bra_last);

        return;
    }

    if (angmom == 1)
    {
        quadrec::compQuadrupolePP(matrix_xx, matrix_xy, matrix_xz, matrix_yy, matrix_yz, matrix_zz, point, gto_block, bra_first, bra_last);

        return;
    }

    if (angmom == 2)
    {
        quadrec::compQuadrupoleDD(matrix_xx, matrix_xy, matrix_xz, matrix_yy, matrix_yz, matrix_zz, point, gto_block, bra_first, bra_last);

        return;
    }

    if (angmom == 3)
    {
        quadrec::compQuadrupoleFF(matrix_xx, matrix_xy, matrix_xz, matrix_yy, matrix_yz, matrix_zz, point, gto_block, bra_first, bra_last);

        return;
    }
}

auto
compute(CSubMatrix*      matrix_xx,
        CSubMatrix*      matrix_xy,
        CSubMatrix*      matrix_xz,
        CSubMatrix*      matrix_yy,
        CSubMatrix*      matrix_yz,
        CSubMatrix*      matrix_zz,
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
        quadrec::compQuadrupoleSS(matrix_xx, matrix_xy, matrix_xz, matrix_yy, matrix_yz, matrix_zz, point, bra_gto_block, ket_gto_block, bra_first, bra_last, mat_type);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 1))
    {
        quadrec::compQuadrupoleSP(matrix_xx, matrix_xy, matrix_xz, matrix_yy, matrix_yz, matrix_zz, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 0))
    {
        quadrec::compQuadrupolePS(matrix_xx, matrix_xy, matrix_xz, matrix_yy, matrix_yz, matrix_zz, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 2))
    {
        quadrec::compQuadrupoleSD(matrix_xx, matrix_xy, matrix_xz, matrix_yy, matrix_yz, matrix_zz, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 0))
    {
        quadrec::compQuadrupoleDS(matrix_xx, matrix_xy, matrix_xz, matrix_yy, matrix_yz, matrix_zz, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 3))
    {
        quadrec::compQuadrupoleSF(matrix_xx, matrix_xy, matrix_xz, matrix_yy, matrix_yz, matrix_zz, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 0))
    {
        quadrec::compQuadrupoleFS(matrix_xx, matrix_xy, matrix_xz, matrix_yy, matrix_yz, matrix_zz, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 1))
    {
        quadrec::compQuadrupolePP(matrix_xx, matrix_xy, matrix_xz, matrix_yy, matrix_yz, matrix_zz, point, bra_gto_block, ket_gto_block, bra_first, bra_last, mat_type);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 2))
    {
        quadrec::compQuadrupolePD(matrix_xx, matrix_xy, matrix_xz, matrix_yy, matrix_yz, matrix_zz, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 1))
    {
        quadrec::compQuadrupoleDP(matrix_xx, matrix_xy, matrix_xz, matrix_yy, matrix_yz, matrix_zz, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 3))
    {
        quadrec::compQuadrupolePF(matrix_xx, matrix_xy, matrix_xz, matrix_yy, matrix_yz, matrix_zz, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 1))
    {
        quadrec::compQuadrupoleFP(matrix_xx, matrix_xy, matrix_xz, matrix_yy, matrix_yz, matrix_zz, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 2))
    {
        quadrec::compQuadrupoleDD(matrix_xx, matrix_xy, matrix_xz, matrix_yy, matrix_yz, matrix_zz, point, bra_gto_block, ket_gto_block, bra_first, bra_last, mat_type);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 3))
    {
        quadrec::compQuadrupoleDF(matrix_xx, matrix_xy, matrix_xz, matrix_yy, matrix_yz, matrix_zz, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 2))
    {
        quadrec::compQuadrupoleFD(matrix_xx, matrix_xy, matrix_xz, matrix_yy, matrix_yz, matrix_zz, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 3))
    {
        quadrec::compQuadrupoleFF(matrix_xx, matrix_xy, matrix_xz, matrix_yy, matrix_yz, matrix_zz, point, bra_gto_block, ket_gto_block, bra_first, bra_last, mat_type);

        return;
    }
}

}  // namespace quadfunc
