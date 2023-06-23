#include "KineticEnergyFunc.hpp"

#include "KineticEnergyRecDD.hpp"
#include "KineticEnergyRecDF.hpp"
#include "KineticEnergyRecDG.hpp"
#include "KineticEnergyRecDP.hpp"
#include "KineticEnergyRecDS.hpp"
#include "KineticEnergyRecFD.hpp"
#include "KineticEnergyRecFF.hpp"
#include "KineticEnergyRecFG.hpp"
#include "KineticEnergyRecFP.hpp"
#include "KineticEnergyRecFS.hpp"
#include "KineticEnergyRecGD.hpp"
#include "KineticEnergyRecGF.hpp"
#include "KineticEnergyRecGG.hpp"
#include "KineticEnergyRecGP.hpp"
#include "KineticEnergyRecGS.hpp"
#include "KineticEnergyRecPD.hpp"
#include "KineticEnergyRecPF.hpp"
#include "KineticEnergyRecPG.hpp"
#include "KineticEnergyRecPP.hpp"
#include "KineticEnergyRecPS.hpp"
#include "KineticEnergyRecSD.hpp"
#include "KineticEnergyRecSF.hpp"
#include "KineticEnergyRecSG.hpp"
#include "KineticEnergyRecSP.hpp"
#include "KineticEnergyRecSS.hpp"


namespace kinfunc {  // kinfunc namespace

auto
compute(CSubMatrix* matrix, const CGtoBlock& gto_block, const int64_t angmom, const int64_t bra_first, const int64_t bra_last) -> void
{
    if (angmom == 0)
    {
        kinrec::compKineticEnergySS(matrix, gto_block, bra_first, bra_last);
        
        return;
    }
    
    if (angmom == 1)
    {
        kinrec::compKineticEnergyPP(matrix, gto_block, bra_first, bra_last);
        
        return;
    }
    
    if (angmom == 2)
    {
        kinrec::compKineticEnergyDD(matrix, gto_block, bra_first, bra_last);
        
        return;
    }
    
    if (angmom == 3)
    {
        kinrec::compKineticEnergyFF(matrix, gto_block, bra_first, bra_last);
        
        return;
    }
    
    if (angmom == 4)
    {
        kinrec::compKineticEnergyGG(matrix, gto_block, bra_first, bra_last);
        
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
        kinrec::compKineticEnergySS(matrix, bra_gto_block, ket_gto_block, bra_first, bra_last, mat_type);
        
        return;
    }
    
    if ((bra_angmom == 0) && (ket_angmom == 1))
    {
        kinrec::compKineticEnergySP(matrix, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);
        
        return;
    }
    
    if ((bra_angmom == 1) && (ket_angmom == 0))
    {
        kinrec::compKineticEnergyPS(matrix, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);
        
        return;
    }
    
    if ((bra_angmom == 0) && (ket_angmom == 2))
    {
        kinrec::compKineticEnergySD(matrix, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);
        
        return;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 0))
    {
        kinrec::compKineticEnergyDS(matrix, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);
        
        return;
    }
    
    if ((bra_angmom == 0) && (ket_angmom == 3))
    {
        kinrec::compKineticEnergySF(matrix, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);
        
        return;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 0))
    {
        kinrec::compKineticEnergyFS(matrix, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);
        
        return;
    }
    
    if ((bra_angmom == 0) && (ket_angmom == 4))
    {
        kinrec::compKineticEnergySG(matrix, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);
        
        return;
    }
    
    if ((bra_angmom == 4) && (ket_angmom == 0))
    {
        kinrec::compKineticEnergyGS(matrix, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);
        
        return;
    }
    
    if ((bra_angmom == 1) && (ket_angmom == 1))
    {
        kinrec::compKineticEnergyPP(matrix, bra_gto_block, ket_gto_block, bra_first, bra_last, mat_type);
        
        return;
    }
    
    if ((bra_angmom == 1) && (ket_angmom == 2))
    {
        kinrec::compKineticEnergyPD(matrix, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);
        
        return;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 1))
    {
        kinrec::compKineticEnergyDP(matrix, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);
        
        return;
    }
    
    if ((bra_angmom == 1) && (ket_angmom == 3))
    {
        kinrec::compKineticEnergyPF(matrix, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);
        
        return;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 1))
    {
        kinrec::compKineticEnergyFP(matrix, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);
        
        return;
    }
    
    if ((bra_angmom == 1) && (ket_angmom == 4))
    {
        kinrec::compKineticEnergyPG(matrix, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);
        
        return;
    }
    
    if ((bra_angmom == 4) && (ket_angmom == 1))
    {
        kinrec::compKineticEnergyGP(matrix, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);
        
        return;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 2))
    {
        kinrec::compKineticEnergyDD(matrix, bra_gto_block, ket_gto_block, bra_first, bra_last, mat_type);
        
        return;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 3))
    {
        kinrec::compKineticEnergyDF(matrix, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);
        
        return;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 2))
    {
        kinrec::compKineticEnergyFD(matrix, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);
        
        return;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 4))
    {
        kinrec::compKineticEnergyDG(matrix, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);
        
        return;
    }
    
    if ((bra_angmom == 4) && (ket_angmom == 2))
    {
        kinrec::compKineticEnergyGD(matrix, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);
        
        return;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 3))
    {
        kinrec::compKineticEnergyFF(matrix, bra_gto_block, ket_gto_block, bra_first, bra_last, mat_type);
        
        return;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 4))
    {
        kinrec::compKineticEnergyFG(matrix, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);
        
        return;
    }
    
    if ((bra_angmom == 4) && (ket_angmom == 3))
    {
        kinrec::compKineticEnergyGF(matrix, bra_gto_block, ket_gto_block, ang_order, bra_first, bra_last);
        
        return;
    }
    
    if ((bra_angmom == 4) && (ket_angmom == 4))
    {
        kinrec::compKineticEnergyGG(matrix, bra_gto_block, ket_gto_block, bra_first, bra_last, mat_type);
        
        return;
    }
}

}  // namespace kinfunc
