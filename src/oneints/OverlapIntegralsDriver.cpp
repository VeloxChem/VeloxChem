//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "OverlapIntegralsDriver.hpp"

#include "OneIntsFunc.hpp"
#include "AngularMomentum.hpp"
#include "OverlapRecFunc.hpp"

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
    
    CMemBlock2D<double> rab(pdim, 3);
    
    auto pmax = bragtos.getMaxContractionDepth();
    
    CMemBlock2D<double> rfacts(pdim, 2 * pmax);
    
    CMemBlock2D<double> rpa(pdim, 3 * pmax);
    
    CMemBlock2D<double> rpb(pdim, 3 * pmax);
    
    // allocate primitives buffer
    
    auto nblk = _getNumberOfPrimBlocks(bragtos, ketgtos, pmax);
    
    CMemBlock2D<double> pbuffer(pdim, nblk);
    
    printf("buffer: %i %i => %i x %i\n", bragtos.getAngularMomentum(),
           ketgtos.getAngularMomentum(), nblk, pdim);

    for (int32_t i = 0; i < bragtos.getNumberOfContrGtos(); i++)
    {
        // compute distances: R(AB) = A - B
        
        intsfunc::compDistancesAB(rab, bragtos, ketgtos, i);
        
        // compute Obara-Saika recursion factors
        
        intsfunc::compFactorsForOverlap(rfacts, bragtos, ketgtos, i);
        
        // compute distances: R(PA) = P - A
        
        intsfunc::compDistancesPA(rpa, rab, rfacts, 2, bragtos, ketgtos, i);
        
        // compute distances: R(PB) = P - B
        
        intsfunc::compDistancesPB(rpb, rab, rfacts, 2, bragtos, ketgtos, i);
        
        // compite primitive overlap integrals
        
        _compPrimOverlapInts(pbuffer, rfacts, rab, rpa, rpb, bragtos, ketgtos, i);
    }
    
    printf("(%i,%i) pair: bra %i ket %i prim (%i,%i) \n", iBraGtoBlock, iKetGtoBlock,
           bragtos.getAngularMomentum(), ketgtos.getAngularMomentum(),
           bragtos.getNumberOfPrimGtos(), ketgtos.getNumberOfPrimGtos());
}

int32_t
COverlapIntegralsDriver::_getNumberOfPrimBlocks(const CGtoBlock& braGtoBlock,
                                                const CGtoBlock& ketGtoBlock,
                                                const int32_t    maxPrimGtos) const
{
    // set up angular momentum data
    
    auto bang = braGtoBlock.getAngularMomentum();
    
    auto kang = ketGtoBlock.getAngularMomentum();
    
    // compute number of blocks in recursion buffer
    
    int32_t ndim = 0;
    
    if (bang > kang)
    {
        // add ket side recursion blocks
        
        for (int32_t i = 0; i <= kang; i++)
        {
            auto kdim = angmom::to_CartesianComponents(i);
            
            for (int32_t j = i; j <= bang; j++)
            {
                ndim += kdim * angmom::to_CartesianComponents(j);
                
                printf("(%i,%i) %i\n", j, i, ndim);
            }
        }
    }
    else
    {
        for (int32_t i = 0; i <= bang; i++)
        {
            auto bdim = angmom::to_CartesianComponents(i);
            
            for (int32_t j = i; j <= kang; j++)
            {
                ndim += bdim * angmom::to_CartesianComponents(j);
                
                printf("(%i,%i) %i\n", i, j, ndim);
            }
        }
    }
    
    return ndim * maxPrimGtos;
}

void
COverlapIntegralsDriver::_compPrimOverlapInts(      CMemBlock2D<double>& primBuffer,
                                              const CMemBlock2D<double>& osFactors,
                                              const CMemBlock2D<double>& abDistances,
                                              const CMemBlock2D<double>& paDistances,
                                              const CMemBlock2D<double>& pbDistances,
                                              const CGtoBlock&           braGtoBlock,
                                              const CGtoBlock&           ketGtoBlock,
                                              const int32_t              iContrGto) const
{
    // set up angular momentum information
    
    auto bang = braGtoBlock.getAngularMomentum();
    
    auto kang = ketGtoBlock.getAngularMomentum();
    
    // compute (s|s) integrals
    
    ovlrecfunc::compOverlapForSS(primBuffer, osFactors, abDistances, braGtoBlock,
                                 ketGtoBlock, iContrGto);
    
    if ((bang == 0) && (kang == 0)) return;
    
    // split Obara-Saika recursion into two alternative paths
    
    if (bang > kang)
    {
        ovlrecfunc::compOverlapForPS(primBuffer, paDistances, braGtoBlock,
                                     ketGtoBlock, iContrGto);
        
        if ((bang == 1) && (kang == 0)) return;
    }
    else
    {
        ovlrecfunc::compOverlapForSP(primBuffer, pbDistances, braGtoBlock,
                                     ketGtoBlock, iContrGto);
        
        if ((bang == 0) && (kang == 1)) return;
        
        ovlrecfunc::compOverlapForSD(primBuffer, osFactors, pbDistances,
                                     braGtoBlock, ketGtoBlock, iContrGto);
        
        if ((bang == 0) && (kang == 2)) return;
        
        ovlrecfunc::compOverlapForSF(primBuffer, osFactors, pbDistances,
                                     braGtoBlock, ketGtoBlock, iContrGto);
        
        if ((bang == 0) && (kang == 3)) return;
        
        ovlrecfunc::compOverlapForSG(primBuffer, osFactors, pbDistances,
                                     braGtoBlock, ketGtoBlock, iContrGto);
        
        if ((bang == 0) && (kang == 4)) return;
    }
}
