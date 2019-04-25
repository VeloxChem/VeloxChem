//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "OverlapIntegralsDriver.hpp"

#include "OneIntsFunc.hpp"
#include "AngularMomentum.hpp"
#include "GenFunc.hpp"
#include "MemBlock.hpp"
#include "StringFormat.hpp"

#include "OverlapRecFuncForSX.hpp"
#include "OverlapRecFuncForPX.hpp"
#include "OverlapRecFuncForDX.hpp"
#include "OverlapRecFuncForFF.hpp"
#include "OverlapRecFuncForFG.hpp"
#include "OverlapRecFuncForGF.hpp"
#include "OverlapRecFuncForGG.hpp"

COverlapIntegralsDriver::COverlapIntegralsDriver(MPI_Comm comm)
{
    _locRank  = mpi::rank(comm);
    
    _locNodes = mpi::nodes(comm);

    mpi::duplicate(comm, &_locComm);
}

COverlapIntegralsDriver::~COverlapIntegralsDriver()
{
    mpi::destroy(&_locComm);
}

COverlapMatrix
COverlapIntegralsDriver::compute(const CMolecule&       molecule,
                                 const CMolecularBasis& basis) const
{
    COverlapMatrix ovlmat;
    
    if (_locRank == mpi::master())
    {
        // set up GTOs container
        
        CGtoContainer bracontr(molecule, basis);
        
        // compute overlap integrals
        
        ovlmat = _compOverlapIntegrals(&bracontr, &bracontr);
    }
    
    return ovlmat;
}

COverlapMatrix
COverlapIntegralsDriver::compute(const CMolecule&       molecule,
                                 const CMolecularBasis& braBasis,
                                 const CMolecularBasis& ketBasis) const
{
    COverlapMatrix ovlmat;
    
    if (_locRank == mpi::master())
    {
        // set up GTOs containers
        
        CGtoContainer bracontr(molecule, braBasis);
        
        CGtoContainer ketcontr(molecule, ketBasis);
        
        // compute overlap integrals
        
        ovlmat = _compOverlapIntegrals(&bracontr, &ketcontr);
    }
    
    return ovlmat;
}

COverlapMatrix
COverlapIntegralsDriver::compute(const CMolecule&       braMolecule,
                                 const CMolecule&       ketMolecule,
                                 const CMolecularBasis& basis) const
{
    COverlapMatrix ovlmat;
    
    if (_locRank == mpi::master())
    {
        // set up GTOs containers
        
        CGtoContainer bracontr(braMolecule, basis);
        
        CGtoContainer ketcontr(ketMolecule, basis);
        
        // compute overlap integrals
        
        ovlmat = _compOverlapIntegrals(&bracontr, &ketcontr);
    }
    
    return ovlmat;
}

COverlapMatrix
COverlapIntegralsDriver::compute(const CMolecule&       braMolecule,
                                 const CMolecule&       ketMolecule,
                                 const CMolecularBasis& braBasis,
                                 const CMolecularBasis& ketBasis) const
{
    COverlapMatrix ovlmat;
    
    if (_locRank == mpi::master())
    {
        // set up GTOs containers
        
        CGtoContainer bracontr(braMolecule, braBasis);
        
        CGtoContainer ketcontr(ketMolecule, ketBasis);
        
        // compute overlap integrals
        
        ovlmat = _compOverlapIntegrals(&bracontr, &ketcontr);
    }
    
    return ovlmat;
}

void
COverlapIntegralsDriver::compute(      double*    intsBatch,
                                 const CGtoBlock& braGtoBlock,
                                 const CGtoBlock& ketGtoBlock) const
{
    // determine dimensions of integrals batch
    
    auto nrow = angmom::to_SphericalComponents(braGtoBlock.getAngularMomentum())
    
              * braGtoBlock.getNumberOfContrGtos();
    
    auto ncol = angmom::to_SphericalComponents(ketGtoBlock.getAngularMomentum())
    
              * ketGtoBlock.getNumberOfContrGtos();
    
    // set up distribution pattern
    
    COneIntsDistribution dist(intsBatch, nrow, ncol, dist1e::batch);
    
    // compute overlap integrals
    
    _compOverlapForGtoBlocks(&dist, braGtoBlock, ketGtoBlock);
}

COverlapMatrix
COverlapIntegralsDriver::_compOverlapIntegrals(const CGtoContainer* braGtoContainer,
                                               const CGtoContainer* ketGtoContainer) const
{
    // check if GTOs containers are same on bra and ket sides
    
    auto symbk = ((*braGtoContainer) == (*ketGtoContainer));
    
    // determine dimensions of overlap matrix
    
    auto nrow = braGtoContainer->getNumberOfAtomicOrbitals();
    
    auto ncol = ketGtoContainer->getNumberOfAtomicOrbitals();
    
    // allocate dense matrix for overlap integrals
    
    CDenseMatrix ovlmat(nrow, ncol);
    
    // set up distributio pattern
    
    dist1e dstyp = (symbk) ? dist1e::symsq : dist1e::rect;
    
    COneIntsDistribution* distpat = new COneIntsDistribution(ovlmat.values(),
                                                             nrow, ncol, dstyp);
    
    // compute overlap integral blocks
    
    #pragma omp parallel shared(braGtoContainer, ketGtoContainer, distpat, symbk)
    {
        #pragma omp single nowait
        {
            // determine number of GTOs blocks in bra/ket sides
            
            auto nbra = braGtoContainer->getNumberOfGtoBlocks();
            
            auto nket = ketGtoContainer->getNumberOfGtoBlocks();
            
            // loop over pairs of GTOs blocks
            
            for (int32_t i = 0; i < nbra; i++)
            {
                auto bgtos = braGtoContainer->getGtoBlock(i);
                
                auto joff = (symbk) ? i : 0;
                
                for (int32_t j = joff; j < nket; j++)
                {
                    #pragma omp task firstprivate(j)
                    {
                        auto kgtos = ketGtoContainer->getGtoBlock(j);
                        
                        _compOverlapForGtoBlocks(distpat, bgtos, kgtos);
                    }
                }
            }
        }
    }
    
    // deallocate distribution pattern
    
    delete distpat;
    
    return COverlapMatrix(ovlmat);
}

void
COverlapIntegralsDriver::_compOverlapForGtoBlocks(      COneIntsDistribution* distPattern,
                                                  const CGtoBlock&            braGtoBlock,
                                                  const CGtoBlock&            ketGtoBlock) const
{
    // copy GTOs blocks for bra and ket sides
    
    auto bragtos = braGtoBlock;
    
    auto ketgtos = ketGtoBlock;
    
    // copy distribution pattern
    
    auto distpat = *distPattern;
    
    // set up angular momentum data
    
    auto bang = bragtos.getAngularMomentum();
    
    auto kang = ketgtos.getAngularMomentum();
    
    // set up spherical angular momentum for bra and ket sides
    
    CSphericalMomentum bmom(bang);
    
    CSphericalMomentum kmom(kang);
    
    // allocate prefactors used in Obara-Saika recursion
    
    auto pdim = ketgtos.getNumberOfPrimGtos();
    
    CMemBlock2D<double> rab(pdim, 3);
    
    auto pmax = bragtos.getMaxContractionDepth();
    
    CMemBlock2D<double> rfacts(pdim, 2 * pmax);
    
    // set up tensors of PA and PB distances
    
    auto btcomps = intsfunc::getNumberOfComponentsInDistancesTensor(bang);
    
    auto ktcomps = intsfunc::getNumberOfComponentsInDistancesTensor(kang);
    
    bool userpa2b = btcomps * ktcomps > 0;
    
    auto rpa = (btcomps > 0) ? CMemBlock2D<double>(pdim, btcomps * pmax) : CMemBlock2D<double>();
    
    auto rpb = (ktcomps > 0) ? CMemBlock2D<double>(pdim, ktcomps * pmax) : CMemBlock2D<double>();
    
    auto rpa2b = (userpa2b) ? CMemBlock2D<double>(pdim, btcomps * ktcomps * pmax) : CMemBlock2D<double>();

    // allocate primitives and auxilary integrals buffer
    
    auto ncart = angmom::to_CartesianComponents(bang, kang);
    
    CMemBlock2D<double> primbuffer(pdim, ncart * pmax);
    
    CMemBlock2D<double> auxbuffer(pdim, pmax);
    
    // allocate contracted Cartesian integrals buffer
    
    auto kdim = ketgtos.getNumberOfContrGtos();
    
    CMemBlock2D<double> cartbuffer(kdim, ncart);
    
    // allocate contracted spherical integrals buffer
    
    auto nspher = angmom::to_SphericalComponents(bang, kang);
    
    CMemBlock2D<double> spherbuffer(kdim, nspher);
    
    // determine bra and ket sides symmetry
    
    bool symbk = (bragtos == ketgtos);

    for (int32_t i = 0; i < bragtos.getNumberOfContrGtos(); i++)
    {
        // compute distances: R(AB) = A - B
        
        intsfunc::compDistancesAB(rab, bragtos, ketgtos, i);
        
        // compute Obara-Saika recursion factors
        
        intsfunc::compFactorsForOverlap(rfacts, bragtos, ketgtos, i);
        
        // compute tensors of distances: R(PA) = P - A
        
        intsfunc::compTensorsPA(rpa, rab, rfacts, 2, bragtos, ketgtos, i);
        
        // compute tensors of distances: R(PB) = P - B
        
        intsfunc::compTensorsPB(rpb, rab, rfacts, 2, bragtos, ketgtos, i);
        
        // compute tensor products: R(PA) x P(PB)
        
        intsfunc::compTensorsProduct(rpa2b, rpa, rpb);
        
        // compite primitive overlap integrals
        
        _compPrimOverlapInts(primbuffer, auxbuffer, rfacts, rab, rpa, rpb, rpa2b, 
                             bragtos, ketgtos, i);
        
        // contract primitive overlap integrals
        
        genfunc::contract(cartbuffer, primbuffer, 0, bragtos, ketgtos, i);
        
        // transform Cartesian to spherical integrals
        
        genfunc::transform(spherbuffer, cartbuffer, bmom, kmom, 0, 0, kdim);
        
        // add batch of integrals to integrals matrix
        
        distpat.distribute(spherbuffer, bragtos, ketgtos, symbk, i);
    }
}

void
COverlapIntegralsDriver::_compPrimOverlapInts(      CMemBlock2D<double>&  primBuffer,
                                                    CMemBlock2D<double>&  auxBuffer,
                                              const CMemBlock2D<double>&  osFactors,
                                              const CMemBlock2D<double>&  abDistances,
                                              const CMemBlock2D<double>&  paDistances,
                                              const CMemBlock2D<double>&  pbDistances,
                                              const CMemBlock2D<double>&  pa2bDistances,
                                              const CGtoBlock&            braGtoBlock,
                                              const CGtoBlock&            ketGtoBlock,
                                              const int32_t               iContrGto) const
{
    // set up angular momentum on bra and ket sides
    
    auto bang = braGtoBlock.getAngularMomentum();
    
    auto kang = ketGtoBlock.getAngularMomentum();
    
    // compute (s||s) auxilary integrals
    
    ovlrecfunc::compOverlapForSS(primBuffer, auxBuffer, osFactors, abDistances,
                                 braGtoBlock, ketGtoBlock, iContrGto);
    
    // compute (s||p) auxilary integrals
    
    if ((bang == 0) && (kang == 1))
    {
        ovlrecfunc::compOverlapForSP(primBuffer, auxBuffer, pbDistances,
                                     braGtoBlock, ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (s||d) auxilary integrals
    
    if ((bang == 0) && (kang == 2))
    {
        ovlrecfunc::compOverlapForSD(primBuffer, auxBuffer, osFactors,
                                     pbDistances, braGtoBlock, ketGtoBlock,
                                     iContrGto);
        
        return;
    }
    
    // compute (s||f) auxilary integrals
    
    if ((bang == 0) && (kang == 3))
    {
        ovlrecfunc::compOverlapForSF(primBuffer, auxBuffer, osFactors,
                                     pbDistances, braGtoBlock, ketGtoBlock,
                                     iContrGto);
        
        return;
    }
    
    // compute (s||g) auxilary integrals
    
    if ((bang == 0) && (kang == 4))
    {
        ovlrecfunc::compOverlapForSG(primBuffer, auxBuffer, osFactors,
                                     pbDistances, braGtoBlock, ketGtoBlock,
                                     iContrGto);
        
        return;
    }
    
    // compute (p||s) auxilary integrals
    
    if ((bang == 1) && (kang == 0))
    {
        ovlrecfunc::compOverlapForPS(primBuffer, auxBuffer, paDistances,
                                     braGtoBlock, ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (d||s) auxilary integrals
    
    if ((bang == 2) && (kang == 0))
    {
        ovlrecfunc::compOverlapForDS(primBuffer, auxBuffer, osFactors,
                                     paDistances, braGtoBlock, ketGtoBlock,
                                     iContrGto);
        
        return;
    }
    
    // compute (f||s) auxilary integrals
    
    if ((bang == 3) && (kang == 0))
    {
        ovlrecfunc::compOverlapForFS(primBuffer, auxBuffer, osFactors,
                                     paDistances, braGtoBlock, ketGtoBlock,
                                     iContrGto);
        
        return;
    }
    
    // compute (g||s) auxilary integrals
    
    if ((bang == 4) && (kang == 0))
    {
        ovlrecfunc::compOverlapForGS(primBuffer, auxBuffer, osFactors,
                                     paDistances, braGtoBlock, ketGtoBlock,
                                     iContrGto);
        
        return;
    }
    
    // compute (p||p) auxilary integrals
    
    if ((bang == 1) && (kang == 1))
    {
        ovlrecfunc::compOverlapForPP(primBuffer, auxBuffer, osFactors,
                                     paDistances, pbDistances, braGtoBlock,
                                     ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (p||d) auxilary integrals
    
    if ((bang == 1) && (kang == 2))
    {
        ovlrecfunc::compOverlapForPD(primBuffer, auxBuffer, osFactors,
                                     paDistances, pbDistances, braGtoBlock,
                                     ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (d||p) auxilary integrals
    
    if ((bang == 2) && (kang == 1))
    {
        ovlrecfunc::compOverlapForDP(primBuffer, auxBuffer, osFactors,
                                     paDistances, pbDistances, braGtoBlock,
                                     ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (p||f) auxilary integrals
    
    if ((bang == 1) && (kang == 3))
    {
        ovlrecfunc::compOverlapForPF(primBuffer, auxBuffer, osFactors,
                                     paDistances, pbDistances, braGtoBlock,
                                     ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (f||p) auxilary integrals
    
    if ((bang == 3) && (kang == 1))
    {
        ovlrecfunc::compOverlapForFP(primBuffer, auxBuffer, osFactors,
                                     paDistances, pbDistances, braGtoBlock,
                                     ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (p||g) auxilary integrals
    
    if ((bang == 1) && (kang == 4))
    {
        ovlrecfunc::compOverlapForPG(primBuffer, auxBuffer, osFactors,
                                     paDistances, pbDistances, braGtoBlock,
                                     ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (g||p) auxilary integrals
    
    if ((bang == 4) && (kang == 1))
    {
        ovlrecfunc::compOverlapForGP(primBuffer, auxBuffer, osFactors,
                                     paDistances, pbDistances, braGtoBlock,
                                     ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (d||d) auxilary integrals
    
    if ((bang == 2) && (kang == 2))
    {
        ovlrecfunc::compOverlapForDD(primBuffer, auxBuffer, osFactors,
                                     paDistances, pbDistances, braGtoBlock,
                                     ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (d||f) auxilary integrals
    
    if ((bang == 2) && (kang == 3))
    {
        ovlrecfunc::compOverlapForDF(primBuffer, auxBuffer, osFactors,
                                     paDistances, pbDistances, braGtoBlock,
                                     ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (f||d) auxilary integrals
    
    if ((bang == 3) && (kang == 2))
    {
        ovlrecfunc::compOverlapForFD(primBuffer, auxBuffer, osFactors,
                                     paDistances, pbDistances, braGtoBlock,
                                     ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (d||g) auxilary integrals
    
    if ((bang == 2) && (kang == 4))
    {
        ovlrecfunc::compOverlapForDG(primBuffer, auxBuffer, osFactors,
                                     paDistances, pbDistances, braGtoBlock,
                                     ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (g||d) auxilary integrals
    
    if ((bang == 4) && (kang == 2))
    {
        ovlrecfunc::compOverlapForGD(primBuffer, auxBuffer, osFactors,
                                     paDistances, pbDistances, braGtoBlock,
                                     ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (f||f) auxilary integrals
    
    if ((bang == 3) && (kang == 3))
    {
        ovlrecfunc::compOverlapForFF(primBuffer, auxBuffer, osFactors,
                                     paDistances, pbDistances, braGtoBlock,
                                     ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (f||g) auxilary integrals
    
    if ((bang == 3) && (kang == 4))
    {
        ovlrecfunc::compOverlapForFG(primBuffer, auxBuffer, osFactors,
                                     paDistances, pbDistances, braGtoBlock,
                                     ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (g||f) auxilary integrals
    
    if ((bang == 4) && (kang == 3))
    {
        ovlrecfunc::compOverlapForGF(primBuffer, auxBuffer, osFactors,
                                     paDistances, pbDistances, braGtoBlock,
                                     ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (g||g) auxilary integrals
    
    if ((bang == 4) && (kang == 4))
    {
        ovlrecfunc::compOverlapForGG(primBuffer, auxBuffer, osFactors,
                                     paDistances, pbDistances, braGtoBlock,
                                     ketGtoBlock, iContrGto);
        
        return;
    }
    
    // NOTE: Add l > 4 terms if needed
}
