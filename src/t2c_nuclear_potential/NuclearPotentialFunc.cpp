#include "NuclearPotentialFunc.hpp"

#include "NuclearPotentialRecDD.hpp"
#include "NuclearPotentialRecDF.hpp"
#include "NuclearPotentialRecDP.hpp"
#include "NuclearPotentialRecDS.hpp"
#include "NuclearPotentialRecFD.hpp"
#include "NuclearPotentialRecFF.hpp"
#include "NuclearPotentialRecFP.hpp"
#include "NuclearPotentialRecFS.hpp"
#include "NuclearPotentialRecPD.hpp"
#include "NuclearPotentialRecPF.hpp"
#include "NuclearPotentialRecPP.hpp"
#include "NuclearPotentialRecPS.hpp"
#include "NuclearPotentialRecSD.hpp"
#include "NuclearPotentialRecSF.hpp"
#include "NuclearPotentialRecSP.hpp"
#include "NuclearPotentialRecSS.hpp"

namespace npotfunc {  // npotfunc namespace

auto
compute(CSubMatrix*      matrix,
        const double     charge,
        const TPoint3D&  point,
        const CGtoBlock& gto_block,
        const int64_t    angmom,
        const int64_t    bra_first,
        const int64_t    bra_last) -> void
{
    if (angmom == 0)
    {
        npotrec::compNuclearPotentialSS(matrix, charge, point, gto_block, bra_first, bra_last);

        return;
    }

    if (angmom == 1)
    {
        npotrec::compNuclearPotentialPP(matrix, charge, point, gto_block, bra_first, bra_last);

        return;
    }

    if (angmom == 2)
    {
        npotrec::compNuclearPotentialDD(matrix, charge, point, gto_block, bra_first, bra_last);

        return;
    }

    if (angmom == 3)
    {
        npotrec::compNuclearPotentialFF(matrix, charge, point, gto_block, bra_first, bra_last);

        return;
    }
}

auto
compute(CSubMatrix*      matrix,
        const double     charge,
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
        npotrec::compNuclearPotentialSS(matrix, charge, point, bra_gto_block, ket_gto_block, bra_first, bra_last, mat_type);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 1))
    {
        npotrec::compNuclearPotentialSP(matrix, charge, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 0))
    {
        npotrec::compNuclearPotentialPS(matrix, charge, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 2))
    {
        npotrec::compNuclearPotentialSD(matrix, charge, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 0))
    {
        npotrec::compNuclearPotentialDS(matrix, charge, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 3))
    {
        npotrec::compNuclearPotentialSF(matrix, charge, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 0))
    {
        npotrec::compNuclearPotentialFS(matrix, charge, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 1))
    {
        npotrec::compNuclearPotentialPP(matrix, charge, point, bra_gto_block, ket_gto_block, bra_first, bra_last, mat_type);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 2))
    {
        npotrec::compNuclearPotentialPD(matrix, charge, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 1))
    {
        npotrec::compNuclearPotentialDP(matrix, charge, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 3))
    {
        npotrec::compNuclearPotentialPF(matrix, charge, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 1))
    {
        npotrec::compNuclearPotentialFP(matrix, charge, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 2))
    {
        npotrec::compNuclearPotentialDD(matrix, charge, point, bra_gto_block, ket_gto_block, bra_first, bra_last, mat_type);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 3))
    {
        npotrec::compNuclearPotentialDF(matrix, charge, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 2))
    {
        npotrec::compNuclearPotentialFD(matrix, charge, point, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 3))
    {
        npotrec::compNuclearPotentialFF(matrix, charge, point, bra_gto_block, ket_gto_block, bra_first, bra_last, mat_type);

        return;
    }
}

}  // namespace npotfunc
