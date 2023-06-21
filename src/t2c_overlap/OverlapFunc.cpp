#include "OverlapFunc.hpp"

#include "OverlapRecSS.hpp"
#include "OverlapRecPP.hpp"
#include "OverlapRecDD.hpp"
#include "OverlapRecFF.hpp"
#include "OverlapRecGG.hpp"
#include "OverlapRecSP.hpp"
#include "OverlapRecPS.hpp"
#include "OverlapRecSD.hpp"
#include "OverlapRecDS.hpp"
#include "OverlapRecSF.hpp"
#include "OverlapRecFS.hpp"
#include "OverlapRecSG.hpp"
#include "OverlapRecGS.hpp"
#include "OverlapRecPD.hpp"
#include "OverlapRecDP.hpp"
#include "OverlapRecPF.hpp"
#include "OverlapRecFP.hpp"
#include "OverlapRecPG.hpp"
#include "OverlapRecGP.hpp"
#include "OverlapRecDF.hpp"
#include "OverlapRecFD.hpp"
#include "OverlapRecDG.hpp"
#include "OverlapRecGD.hpp"
#include "OverlapRecFG.hpp"
#include "OverlapRecGF.hpp"

namespace ovlfunc { // ovlfunc namespace

auto
compute(      CSubMatrix* matrix,
        const CGtoBlock&  gto_block,
        const int64_t     angmom,
        const int64_t     bra_first,
        const int64_t     bra_last) -> void
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
    
    if (angmom == 4)
    {
        ovlrec::compOverlapGG(matrix, gto_block, bra_first, bra_last);
        
        return;
    }
}

auto
compute(      CSubMatrix* matrix,
        const CGtoBlock&  bra_gto_block,
        const CGtoBlock&  ket_gto_block,
        const int64_t     bra_angmom,
        const int64_t     ket_angmom,
        const bool        ang_order,
        const int64_t     bra_first,
        const int64_t     bra_last,
        const mat_t       mat_type) -> void
{
    if ((bra_angmom == 0) && (ket_angmom == 0))
    {
        ovlrec::compOverlapSS(matrix, bra_gto_block, ket_gto_block,
                              bra_first, bra_last, mat_type);
        
        return;
    }
    
    if ((bra_angmom == 0) && (ket_angmom == 1))
    {
        ovlrec::compOverlapSP(matrix, bra_gto_block, ket_gto_block,
                              ang_order, bra_first, bra_last);
        
        return;
    }
    
    if ((bra_angmom == 1) && (ket_angmom == 0))
    {
        ovlrec::compOverlapPS(matrix, bra_gto_block, ket_gto_block,
                              ang_order, bra_first, bra_last);
        
        return;
    }
    
    if ((bra_angmom == 0) && (ket_angmom == 2))
    {
        ovlrec::compOverlapSD(matrix, bra_gto_block, ket_gto_block,
                              ang_order, bra_first, bra_last);
        
        return;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 0))
    {
        ovlrec::compOverlapDS(matrix, bra_gto_block, ket_gto_block,
                              ang_order, bra_first, bra_last);
        
        return;
    }
    
    if ((bra_angmom == 0) && (ket_angmom == 3))
    {
        ovlrec::compOverlapSF(matrix, bra_gto_block, ket_gto_block,
                              ang_order, bra_first, bra_last);
        
        return;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 0))
    {
        ovlrec::compOverlapFS(matrix, bra_gto_block, ket_gto_block,
                              ang_order, bra_first, bra_last);
        
        return;
    }
    
    if ((bra_angmom == 0) && (ket_angmom == 4))
    {
        ovlrec::compOverlapSG(matrix, bra_gto_block, ket_gto_block,
                              ang_order, bra_first, bra_last);
        
        return;
    }
    
    if ((bra_angmom == 4) && (ket_angmom == 0))
    {
        ovlrec::compOverlapGS(matrix, bra_gto_block, ket_gto_block,
                              ang_order, bra_first, bra_last);
        
        return;
    }
    
    if ((bra_angmom == 1) && (ket_angmom == 1))
    {
        ovlrec::compOverlapPP(matrix, bra_gto_block, ket_gto_block,
                              bra_first, bra_last, mat_type);
        
        return;
    }
    
    if ((bra_angmom == 1) && (ket_angmom == 2))
    {
        ovlrec::compOverlapPD(matrix, bra_gto_block, ket_gto_block,
                              ang_order, bra_first, bra_last);
        
        return;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 1))
    {
        ovlrec::compOverlapDP(matrix, bra_gto_block, ket_gto_block,
                              ang_order, bra_first, bra_last);
        
        return;
    }
    
    if ((bra_angmom == 1) && (ket_angmom == 3))
    {
        ovlrec::compOverlapPF(matrix, bra_gto_block, ket_gto_block,
                              ang_order, bra_first, bra_last);
        
        return;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 1))
    {
        ovlrec::compOverlapFP(matrix, bra_gto_block, ket_gto_block,
                              ang_order, bra_first, bra_last);
        
        return;
    }
    
    if ((bra_angmom == 1) && (ket_angmom == 4))
    {
        ovlrec::compOverlapPG(matrix, bra_gto_block, ket_gto_block,
                              ang_order, bra_first, bra_last);
        
        return;
    }
    
    if ((bra_angmom == 4) && (ket_angmom == 1))
    {
        ovlrec::compOverlapGP(matrix, bra_gto_block, ket_gto_block,
                              ang_order, bra_first, bra_last);
        
        return;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 2))
    {
        ovlrec::compOverlapDD(matrix, bra_gto_block, ket_gto_block,
                              bra_first, bra_last, mat_type);
        
        return;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 3))
    {
        ovlrec::compOverlapDF(matrix, bra_gto_block, ket_gto_block,
                              ang_order, bra_first, bra_last);
        
        return;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 2))
    {
        ovlrec::compOverlapFD(matrix, bra_gto_block, ket_gto_block,
                              ang_order, bra_first, bra_last);
        
        return;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 4))
    {
        ovlrec::compOverlapDG(matrix, bra_gto_block, ket_gto_block,
                              ang_order, bra_first, bra_last);
        
        return;
    }
    
    if ((bra_angmom == 4) && (ket_angmom == 2))
    {
        ovlrec::compOverlapGD(matrix, bra_gto_block, ket_gto_block,
                              ang_order, bra_first, bra_last);
        
        return;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 3))
    {
        ovlrec::compOverlapFF(matrix, bra_gto_block, ket_gto_block,
                              bra_first, bra_last, mat_type);
        
        return;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 4))
    {
        ovlrec::compOverlapFG(matrix, bra_gto_block, ket_gto_block,
                              ang_order, bra_first, bra_last);
        
        return;
    }
    
    if ((bra_angmom == 4) && (ket_angmom == 3))
    {
        ovlrec::compOverlapGF(matrix, bra_gto_block, ket_gto_block,
                              ang_order, bra_first, bra_last);
        
        return;
    }
    
    if ((bra_angmom == 4) && (ket_angmom == 4))
    {
        ovlrec::compOverlapGG(matrix, bra_gto_block, ket_gto_block,
                              bra_first, bra_last, mat_type);
        
        return;
    }
}

}  // ovlfunc namespace
