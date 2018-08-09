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
