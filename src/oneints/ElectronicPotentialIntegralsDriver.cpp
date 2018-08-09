//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "ElectronicPotentialIntegralsDriver.hpp"

#include "GenFunc.hpp"

CElectronicPotentialIntegralsDriver::CElectronicPotentialIntegralsDriver(const int32_t  globRank,
                                                                         const int32_t  globNodes,
                                                                               MPI_Comm comm)

    : _globRank(globRank)

    , _globNodes(globNodes)

    , _isLocalMode(false)
{
    _locRank  = mpi::rank(comm);
    
    _locNodes = mpi::nodes(comm);
    
    _isLocalMode = !mpi::compare(comm, MPI_COMM_WORLD);
}

CElectronicPotentialIntegralsDriver::~CElectronicPotentialIntegralsDriver()
{
    
}

CElectronicPotentialMatrix
CElectronicPotentialIntegralsDriver::compute(const CMolecule&       molecule,
                                             const CMolecularBasis& basis,
                                                   MPI_Comm         comm) const
{
    if (_locRank == mpi::master())
    {
        // set up GTOs container
        
        CGtoContainer bracontr(molecule, basis);
        
        // compute electronic potential integrals
        
        return _compElectronicPotentialIntegrals(&bracontr, &bracontr);
    }
    
    return CElectronicPotentialMatrix();
}

CElectronicPotentialMatrix
CElectronicPotentialIntegralsDriver::_compElectronicPotentialIntegrals(const CGtoContainer* braGtoContainer,
                                                                       const CGtoContainer* ketGtoContainer) const
{
    // allocate buffer of sparse matrices
    
    auto matbuff = _createSparseBuffer(braGtoContainer, ketGtoContainer);
    
    // compute kinetic energy integral blocks
    
    #pragma omp parallel shared(braGtoContainer, ketGtoContainer, matbuff)
    {
        #pragma omp single nowait
        {
            // determine number of GTOs blocks in bra/ket sides
            
            auto nbra = braGtoContainer->getNumberOfGtoBlocks();
            
            auto nket = ketGtoContainer->getNumberOfGtoBlocks();
            
            // loop over pairs of GTOs blocks
            
            for (int32_t i = 0; i < nbra; i++)
            {
                for (int32_t j = 0; j < nket; j++)
                {
                    #pragma omp task firstprivate(i, j)
                    {
                        _compElectronicPotentialForGtoBlocks(matbuff, braGtoContainer,
                                                             i, ketGtoContainer, j);
                    }
                }
            }
        }
    }
    
    // distribute submatrices into single sparse matrix
    
    auto spmat = genfunc::distribute(matbuff, braGtoContainer, ketGtoContainer);
    
    // optimize memory usage in sparsse matrix
    
    spmat.optimize_storage();
    
    // deallocate buffer of sparse matrices
    
    delete [] matbuff;
    
    return CElectronicPotentialMatrix(spmat);
}

void
CElectronicPotentialIntegralsDriver::_compElectronicPotentialForGtoBlocks(      CSparseMatrix* sparseBuffer,
                                                                          const CGtoContainer* braGtoContainer,
                                                                          const int32_t        iBraGtoBlock,
                                                                          const CGtoContainer* ketGtoContainer,
                                                                          const int32_t        iKetGtoBlock) const
{
    // copy GTOs blocks for bra and ket sides
    
    auto bragtos = braGtoContainer->getGtoBlock(iBraGtoBlock);
    
    auto ketgtos = ketGtoContainer->getGtoBlock(iKetGtoBlock);
    
    // set up spherical angular momentum for bra and ket sides
    
    CSphericalMomentum bmom(bragtos.getAngularMomentum());
    
    CSphericalMomentum kmom(ketgtos.getAngularMomentum());
    
    // allocate prefactors used in Obara-Saika recursion
    
    auto pdim = ketgtos.getNumberOfPrimGtos();
    
    CMemBlock2D<double> rab(pdim, 3);
    
    auto pmax = bragtos.getMaxContractionDepth();
    
    CMemBlock2D<double> rfacts(pdim, 6 * pmax);
    
    CMemBlock2D<double> rpa(pdim, 3 * pmax);
    
    CMemBlock2D<double> rpb(pdim, 3 * pmax);
    
    // generate recursion pattern
    
    auto recvec = _getRecursionPattern(bragtos, ketgtos); 
}

CSparseMatrix*
CElectronicPotentialIntegralsDriver::_createSparseBuffer(const CGtoContainer* braGtoContainer,
                                                         const CGtoContainer* ketGtoContainer) const
{
    // setup size of sparse matrices buffer
    
    auto bcomp = braGtoContainer->getNumberOfGtoBlocks();
    
    auto kcomp = ketGtoContainer->getNumberOfGtoBlocks();
    
    // allocate sparse matrices
    
    CSparseMatrix* matbuff = new CSparseMatrix[bcomp * kcomp];
    
    return matbuff;
}

CVecThreeIndexes
CElectronicPotentialIntegralsDriver::_getRecursionPattern(const CGtoBlock& braGtoBlock,
                                                          const CGtoBlock& ketGtoBlock) const
{
    // set up angular momentum
    
    auto bang = braGtoBlock.getAngularMomentum();
    
    auto kang = ketGtoBlock.getAngularMomentum();
    
    // set up recursion buffer
    
    CVecThreeIndexes recvec;
    
    recvec.reserve((bang + 1) * (kang + 1));
    
    // set up indexing counters
    
    int32_t spos = 0;
    
    int32_t epos = 1;
    
    // set up initial state of recursion buffer
    
    recvec.push_back(CThreeIndexes(bang, kang, 0));
    
    while (true)
    {
        // internal new recursion terms counter
        
        int32_t nterms = 0;
        
        // generate bra and ket Obara-Saika recursion terms
        
        for (int32_t i = spos; i < epos; i++)
        {
            CThreeIndexes cidx(recvec[i]);
            
            // nuclear potentil recursion
            
            if (cidx.first() != 0)
            {
                // general recursion for bra and ket sides
                
                // (a - 1 |A(0)| b)^(m) term
                
                CThreeIndexes t10idx(cidx.first() - 1,  cidx.second(),
                                     
                                     cidx.third());
                
                if (genfunc::addValidAndUniqueTriple(recvec, t10idx)) nterms++;
                
                // (a - 1 |A(0)| b)^(m+1) term
                
                CThreeIndexes t11idx(cidx.first() - 1,  cidx.second(),
                                     
                                     cidx.third() + 1);
                
                if (genfunc::addValidAndUniqueTriple(recvec, t11idx)) nterms++;
                
                // (a - 2 |A(0)| b)^(m) term
                
                CThreeIndexes t20idx(cidx.first() - 2,  cidx.second(),
                                     
                                     cidx.third());
                
                if (genfunc::addValidAndUniqueTriple(recvec, t20idx)) nterms++;
                
                // (a - 2 |A(0)| b)^(m+1) term
                
                CThreeIndexes t21idx(cidx.first() - 2,  cidx.second(),
                                     
                                     cidx.third() + 1);
                
                if (genfunc::addValidAndUniqueTriple(recvec, t21idx)) nterms++;
                
                // (a - 1 |A(0)| b - 1)^(m) term
                
                CThreeIndexes tk0idx(cidx.first() - 1,  cidx.second() - 1,
                                     
                                     cidx.third());
                
                if (genfunc::addValidAndUniqueTriple(recvec, tk0idx)) nterms++;
                
                // (a - 1 |A(0)| b - 1)^(m+1) term
                
                CThreeIndexes tk1idx(cidx.first() - 1,  cidx.second() - 1,
                                     
                                     cidx.third() + 1);
                
                if (genfunc::addValidAndUniqueTriple(recvec, tk1idx)) nterms++;
            }
            else
            {
                // simplified recursion for ket sides
                
                // (0 |A(0)| b - 1)^(m) term
                
                CThreeIndexes t10idx(cidx.first(),  cidx.second() - 1,
                                     
                                     cidx.third());
                
                if (genfunc::addValidAndUniqueTriple(recvec, t10idx)) nterms++;
                
                // (0 |A(0)| b - 1)^(m+1) term
                
                CThreeIndexes t11idx(cidx.first(),  cidx.second() - 1,
                                     
                                     cidx.third() + 1);
                
                if (genfunc::addValidAndUniqueTriple(recvec, t11idx)) nterms++;
                
                // (0 |A(0)| b - 2)^(m) term
                
                CThreeIndexes t20idx(cidx.first(),  cidx.second() - 2,
                                     
                                     cidx.third());
                
                if (genfunc::addValidAndUniqueTriple(recvec, t20idx)) nterms++;
                
                // (0 |A(0)| b - 2)^(m+1) term
                
                CThreeIndexes t21idx(cidx.first(),  cidx.second() - 2,
                                     
                                     cidx.third() + 1);
                
                if (genfunc::addValidAndUniqueTriple(recvec, t21idx)) nterms++;
            }
        }
        
        // break loop, all recursion terms are generrated
        
        if (nterms == 0) break;
        
        // update counters
        
        spos  = epos;
        
        epos += nterms;
    }
    
    return recvec;
}
