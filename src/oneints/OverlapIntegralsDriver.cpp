//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "OverlapIntegralsDriver.hpp"

#include "OneIntsFunc.hpp"

COverlapIntegralsDriver::COverlapIntegralsDriver(const int32_t  globRank,
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

COverlapIntegralsDriver::~COverlapIntegralsDriver()
{
}

COverlapMatrix
COverlapIntegralsDriver::compute(const CMolecule&       molecule,
                                 const CMolecularBasis& basis,
                                       MPI_Comm         comm) const 
{
    if (_locRank == mpi::master())
    {
        // set up GTOs container
        
        CGtoContainer bracontr(molecule, basis);
        
        // compute overlap integrals
        
        _compOverlapIntegrals(&bracontr, &bracontr);
    }
    
    return COverlapMatrix();
}

COverlapMatrix
COverlapIntegralsDriver::compute(const CMolecule&       molecule,
                                 const CMolecularBasis& braBasis,
                                 const CMolecularBasis& ketBasis,
                                       MPI_Comm         comm) const
{
    if (_locRank == mpi::master())
    {
        // set up GTOs containers
        
        CGtoContainer bracontr(molecule, braBasis);
        
        CGtoContainer ketcontr(molecule, ketBasis);
        
        // compute overlap integrals
        
        _compOverlapIntegrals(&bracontr, &ketcontr);
    }
    
    return COverlapMatrix();
}

COverlapMatrix
COverlapIntegralsDriver::compute(const CMolecule&       braMolecule,
                                 const CMolecule&       ketMolecule,
                                 const CMolecularBasis& basis,
                                       MPI_Comm         comm) const
{
    if (_locRank == mpi::master())
    {
        // set up GTOs containers
        
        CGtoContainer bracontr(braMolecule, basis);
        
        CGtoContainer ketcontr(ketMolecule, basis);
        
        // compute overlap integrals
        
        _compOverlapIntegrals(&bracontr, &ketcontr);
    }
    
    return COverlapMatrix();
}

COverlapMatrix
COverlapIntegralsDriver::compute(const CMolecule&       braMolecule,
                                 const CMolecule&       ketMolecule,
                                 const CMolecularBasis& braBasis,
                                 const CMolecularBasis& ketBasis,
                                       MPI_Comm         comm) const
{
    if (_locRank == mpi::master())
    {
        // set up GTOs containers
        
        CGtoContainer bracontr(braMolecule, braBasis);
        
        CGtoContainer ketcontr(ketMolecule, ketBasis);
        
        // compute overlap integrals
        
        _compOverlapIntegrals(&bracontr, &ketcontr);
    }
    
    return COverlapMatrix();
}

void
COverlapIntegralsDriver::_compOverlapIntegrals(const CGtoContainer* braGtoContainer,
                                               const CGtoContainer* ketGtoContainer) const
{
    // compute overlap integral blocks
    
    #pragma omp parallel shared(braGtoContainer, ketGtoContainer)
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
                        _compOverlapForGtoBlocks(braGtoContainer, i,
                                                 ketGtoContainer, j);
                    }
                }
            }
        }
    }
}


void
COverlapIntegralsDriver::_compOverlapForGtoBlocks(const CGtoContainer* braGtoContainer,
                                                  const int32_t        iBraGtoBlock,
                                                  const CGtoContainer* ketGtoContainer,
                                                  const int32_t        iKetGtoBlock) const
{
    // copy GTOs blocks for bra and ket sides
    
    auto bragtos = braGtoContainer->getGtoBlock(iBraGtoBlock);
    
    auto ketgtos = ketGtoContainer->getGtoBlock(iKetGtoBlock);
    
    // allocate prefactors used in Obara-Saika recursion
    
    auto pdim = ketgtos.getNumberOfPrimGtos();
    
    auto pmax = bragtos.getMaxContractionDepth();
    
    CMemBlock2D<double> rab(pdim, 3);
    
    CMemBlock2D<double> rfacts(pdim, 2 * pmax);
    
    CMemBlock2D<double> rpa(pdim, 3 * pmax);
    
    CMemBlock2D<double> rpb(pdim, 3 * pmax);
    
    // set up pointers to primitives data in bra side
    
    auto brx = bragtos.getCoordinatesX();
    
    auto bry = bragtos.getCoordinatesY();
    
    auto brz = bragtos.getCoordinatesZ();
    
    auto bexp = bragtos.getExponents();
    
    auto bang = bragtos.getAngularMomentum();
    
    // set up pointets to primitives data in ket side
    
    auto krx = ketgtos.getCoordinatesX();
    
    auto kry = ketgtos.getCoordinatesY();
    
    auto krz = ketgtos.getCoordinatesZ();
    
    auto kexp = ketgtos.getExponents();
    
    auto kang = ketgtos.getAngularMomentum();
    
    // loop over contracted functions in bra side
    
    auto spos = bragtos.getStartPositions();
    
    auto epos = bragtos.getEndPositions();
    
    for (int32_t i = 0; i < bragtos.getNumberOfContrGtos(); i++)
    {
        // compute various prefactors for contracted GTO
        
        int32_t idx = 0;
        
        for (int32_t j = spos[i]; j < epos[i]; j++)
        {
            // compute distances: R(A-B)
            
            if (idx == 0)
            {
                mathfunc::distances(rab.data(0), rab.data(1), rab.data(2),
                                    brx[j], bry[j], brz[j], krx, kry, krz,
                                    pdim);
            }
            
            // compute Xi and Zeta factors
            
            intsfunc::compXiAndZeta(rfacts.data(2 * idx), rfacts.data(2 * idx + 1),
                                    bexp[j], kexp, pdim);
            
            // compute P-A distances
            
            if (bang > 0)
            {
                intsfunc::compDistancesPA(rpa.data(3 * idx), rpa.data(3 * idx + 1),
                                          rpa.data(3 * idx + 2), kexp,
                                          rfacts.data(2 * idx), rab.data(0),
                                          rab.data(1), rab.data(2), pdim);
            }
            
            // compute P-B distances
            
            if (kang > 0)
            {
                intsfunc::compDistancesPB(rpb.data(3 * idx), rpb.data(3 * idx + 1),
                                          rpb.data(3 * idx + 2), bexp[j],
                                          rfacts.data(2 * idx), rab.data(0),
                                          rab.data(1), rab.data(2), pdim);
            }
            
            // compute batch of primitive overlap integrals
            
            //_compPrimOverlap()
            
            idx++;
        }
        
        // contract primitive overlap integrals
        
    }
    
    printf("(%i,%i) pair: bra %i ket %i prim (%i,%i) \n", iBraGtoBlock, iKetGtoBlock,
           bragtos.getAngularMomentum(), ketgtos.getAngularMomentum(),
           bragtos.getNumberOfPrimGtos(), ketgtos.getNumberOfPrimGtos());
}
