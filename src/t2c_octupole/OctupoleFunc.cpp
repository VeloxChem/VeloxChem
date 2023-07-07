#include "OctupoleFunc.hpp"

#include "OctupoleRecDD.hpp"
#include "OctupoleRecDF.hpp"
#include "OctupoleRecDP.hpp"
#include "OctupoleRecDS.hpp"
#include "OctupoleRecFD.hpp"
#include "OctupoleRecFF.hpp"
#include "OctupoleRecFP.hpp"
#include "OctupoleRecFS.hpp"
#include "OctupoleRecPD.hpp"
#include "OctupoleRecPF.hpp"
#include "OctupoleRecPP.hpp"
#include "OctupoleRecPS.hpp"
#include "OctupoleRecSD.hpp"
#include "OctupoleRecSF.hpp"
#include "OctupoleRecSP.hpp"
#include "OctupoleRecSS.hpp"

namespace octufunc {  // octufunc namespace

auto
compute(CSubMatrix*      matrix_xxx,
        CSubMatrix*      matrix_xxy,
        CSubMatrix*      matrix_xxz,
        CSubMatrix*      matrix_xyy,
        CSubMatrix*      matrix_xyz,
        CSubMatrix*      matrix_xzz,
        CSubMatrix*      matrix_yyy,
        CSubMatrix*      matrix_yyz,
        CSubMatrix*      matrix_yzz,
        CSubMatrix*      matrix_zzz,
        const TPoint3D&  point,
        const CGtoBlock& gto_block,
        const int64_t    angmom,
        const int64_t    bra_first,
        const int64_t    bra_last) -> void
{
    if (angmom == 0)
    {
        octurec::compOctupoleSS(matrix_xxx,
                                matrix_xxy,
                                matrix_xxz,
                                matrix_xyy,
                                matrix_xyz,
                                matrix_xzz,
                                matrix_yyy,
                                matrix_yyz,
                                matrix_yzz,
                                matrix_zzz,
                                point,
                                gto_block,
                                bra_first,
                                bra_last);

        return;
    }

    if (angmom == 1)
    {
        octurec::compOctupolePP(matrix_xxx,
                                matrix_xxy,
                                matrix_xxz,
                                matrix_xyy,
                                matrix_xyz,
                                matrix_xzz,
                                matrix_yyy,
                                matrix_yyz,
                                matrix_yzz,
                                matrix_zzz,
                                point,
                                gto_block,
                                bra_first,
                                bra_last);

        return;
    }

    if (angmom == 2)
    {
        octurec::compOctupoleDD(matrix_xxx,
                                matrix_xxy,
                                matrix_xxz,
                                matrix_xyy,
                                matrix_xyz,
                                matrix_xzz,
                                matrix_yyy,
                                matrix_yyz,
                                matrix_yzz,
                                matrix_zzz,
                                point,
                                gto_block,
                                bra_first,
                                bra_last);

        return;
    }

    if (angmom == 3)
    {
        octurec::compOctupoleFF(matrix_xxx,
                                matrix_xxy,
                                matrix_xxz,
                                matrix_xyy,
                                matrix_xyz,
                                matrix_xzz,
                                matrix_yyy,
                                matrix_yyz,
                                matrix_yzz,
                                matrix_zzz,
                                point,
                                gto_block,
                                bra_first,
                                bra_last);

        return;
    }
}

auto
compute(CSubMatrix*      matrix_xxx,
        CSubMatrix*      matrix_xxy,
        CSubMatrix*      matrix_xxz,
        CSubMatrix*      matrix_xyy,
        CSubMatrix*      matrix_xyz,
        CSubMatrix*      matrix_xzz,
        CSubMatrix*      matrix_yyy,
        CSubMatrix*      matrix_yyz,
        CSubMatrix*      matrix_yzz,
        CSubMatrix*      matrix_zzz,
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
        octurec::compOctupoleSS(matrix_xxx,
                                matrix_xxy,
                                matrix_xxz,
                                matrix_xyy,
                                matrix_xyz,
                                matrix_xzz,
                                matrix_yyy,
                                matrix_yyz,
                                matrix_yzz,
                                matrix_zzz,
                                point,
                                bra_gto_block,
                                ket_gto_block,
                                bra_first,
                                bra_last,
                                mat_type);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 1))
    {
        octurec::compOctupoleSP(matrix_xxx,
                                matrix_xxy,
                                matrix_xxz,
                                matrix_xyy,
                                matrix_xyz,
                                matrix_xzz,
                                matrix_yyy,
                                matrix_yyz,
                                matrix_yzz,
                                matrix_zzz,
                                point,
                                bra_gto_block,
                                ket_gto_block,
                                ang_order,
                                bra_first,
                                bra_last);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 0))
    {
        octurec::compOctupolePS(matrix_xxx,
                                matrix_xxy,
                                matrix_xxz,
                                matrix_xyy,
                                matrix_xyz,
                                matrix_xzz,
                                matrix_yyy,
                                matrix_yyz,
                                matrix_yzz,
                                matrix_zzz,
                                point,
                                bra_gto_block,
                                ket_gto_block,
                                ang_order,
                                bra_first,
                                bra_last);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 2))
    {
        octurec::compOctupoleSD(matrix_xxx,
                                matrix_xxy,
                                matrix_xxz,
                                matrix_xyy,
                                matrix_xyz,
                                matrix_xzz,
                                matrix_yyy,
                                matrix_yyz,
                                matrix_yzz,
                                matrix_zzz,
                                point,
                                bra_gto_block,
                                ket_gto_block,
                                ang_order,
                                bra_first,
                                bra_last);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 0))
    {
        octurec::compOctupoleDS(matrix_xxx,
                                matrix_xxy,
                                matrix_xxz,
                                matrix_xyy,
                                matrix_xyz,
                                matrix_xzz,
                                matrix_yyy,
                                matrix_yyz,
                                matrix_yzz,
                                matrix_zzz,
                                point,
                                bra_gto_block,
                                ket_gto_block,
                                ang_order,
                                bra_first,
                                bra_last);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 3))
    {
        octurec::compOctupoleSF(matrix_xxx,
                                matrix_xxy,
                                matrix_xxz,
                                matrix_xyy,
                                matrix_xyz,
                                matrix_xzz,
                                matrix_yyy,
                                matrix_yyz,
                                matrix_yzz,
                                matrix_zzz,
                                point,
                                bra_gto_block,
                                ket_gto_block,
                                ang_order,
                                bra_first,
                                bra_last);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 0))
    {
        octurec::compOctupoleFS(matrix_xxx,
                                matrix_xxy,
                                matrix_xxz,
                                matrix_xyy,
                                matrix_xyz,
                                matrix_xzz,
                                matrix_yyy,
                                matrix_yyz,
                                matrix_yzz,
                                matrix_zzz,
                                point,
                                bra_gto_block,
                                ket_gto_block,
                                ang_order,
                                bra_first,
                                bra_last);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 1))
    {
        octurec::compOctupolePP(matrix_xxx,
                                matrix_xxy,
                                matrix_xxz,
                                matrix_xyy,
                                matrix_xyz,
                                matrix_xzz,
                                matrix_yyy,
                                matrix_yyz,
                                matrix_yzz,
                                matrix_zzz,
                                point,
                                bra_gto_block,
                                ket_gto_block,
                                bra_first,
                                bra_last,
                                mat_type);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 2))
    {
        octurec::compOctupolePD(matrix_xxx,
                                matrix_xxy,
                                matrix_xxz,
                                matrix_xyy,
                                matrix_xyz,
                                matrix_xzz,
                                matrix_yyy,
                                matrix_yyz,
                                matrix_yzz,
                                matrix_zzz,
                                point,
                                bra_gto_block,
                                ket_gto_block,
                                ang_order,
                                bra_first,
                                bra_last);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 1))
    {
        octurec::compOctupoleDP(matrix_xxx,
                                matrix_xxy,
                                matrix_xxz,
                                matrix_xyy,
                                matrix_xyz,
                                matrix_xzz,
                                matrix_yyy,
                                matrix_yyz,
                                matrix_yzz,
                                matrix_zzz,
                                point,
                                bra_gto_block,
                                ket_gto_block,
                                ang_order,
                                bra_first,
                                bra_last);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 3))
    {
        octurec::compOctupolePF(matrix_xxx,
                                matrix_xxy,
                                matrix_xxz,
                                matrix_xyy,
                                matrix_xyz,
                                matrix_xzz,
                                matrix_yyy,
                                matrix_yyz,
                                matrix_yzz,
                                matrix_zzz,
                                point,
                                bra_gto_block,
                                ket_gto_block,
                                ang_order,
                                bra_first,
                                bra_last);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 1))
    {
        octurec::compOctupoleFP(matrix_xxx,
                                matrix_xxy,
                                matrix_xxz,
                                matrix_xyy,
                                matrix_xyz,
                                matrix_xzz,
                                matrix_yyy,
                                matrix_yyz,
                                matrix_yzz,
                                matrix_zzz,
                                point,
                                bra_gto_block,
                                ket_gto_block,
                                ang_order,
                                bra_first,
                                bra_last);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 2))
    {
        octurec::compOctupoleDD(matrix_xxx,
                                matrix_xxy,
                                matrix_xxz,
                                matrix_xyy,
                                matrix_xyz,
                                matrix_xzz,
                                matrix_yyy,
                                matrix_yyz,
                                matrix_yzz,
                                matrix_zzz,
                                point,
                                bra_gto_block,
                                ket_gto_block,
                                bra_first,
                                bra_last,
                                mat_type);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 3))
    {
        octurec::compOctupoleDF(matrix_xxx,
                                matrix_xxy,
                                matrix_xxz,
                                matrix_xyy,
                                matrix_xyz,
                                matrix_xzz,
                                matrix_yyy,
                                matrix_yyz,
                                matrix_yzz,
                                matrix_zzz,
                                point,
                                bra_gto_block,
                                ket_gto_block,
                                ang_order,
                                bra_first,
                                bra_last);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 2))
    {
        octurec::compOctupoleFD(matrix_xxx,
                                matrix_xxy,
                                matrix_xxz,
                                matrix_xyy,
                                matrix_xyz,
                                matrix_xzz,
                                matrix_yyy,
                                matrix_yyz,
                                matrix_yzz,
                                matrix_zzz,
                                point,
                                bra_gto_block,
                                ket_gto_block,
                                ang_order,
                                bra_first,
                                bra_last);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 3))
    {
        octurec::compOctupoleFF(matrix_xxx,
                                matrix_xxy,
                                matrix_xxz,
                                matrix_xyy,
                                matrix_xyz,
                                matrix_xzz,
                                matrix_yyy,
                                matrix_yyz,
                                matrix_yzz,
                                matrix_zzz,
                                point,
                                bra_gto_block,
                                ket_gto_block,
                                bra_first,
                                bra_last,
                                mat_type);

        return;
    }
}

}  // namespace octufunc
