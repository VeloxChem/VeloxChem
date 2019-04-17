//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "KineticEnergyIntegralsDriver.hpp"

#include "GenFunc.hpp"
#include "AngularMomentum.hpp"
#include "OneIntsFunc.hpp"
#include "StringFormat.hpp"

#include "KineticEnergyRecFuncForSX.hpp"
#include "KineticEnergyRecFuncForPX.hpp"
#include "KineticEnergyRecFuncForDX.hpp"
#include "KineticEnergyRecFuncForFF.hpp"
#include "KineticEnergyRecFuncForFG.hpp"
#include "KineticEnergyRecFuncForGF.hpp"
#include "KineticEnergyRecFuncForGG.hpp"

CKineticEnergyIntegralsDriver::CKineticEnergyIntegralsDriver(MPI_Comm comm)
{
    _locRank  = mpi::rank(comm);
    
    _locNodes = mpi::nodes(comm);

    mpi::duplicate(comm, &_locComm);
}

CKineticEnergyIntegralsDriver::~CKineticEnergyIntegralsDriver()
{
    mpi::destroy(&_locComm);
}

CKineticEnergyMatrix
CKineticEnergyIntegralsDriver::compute(const CMolecule&       molecule,
                                       const CMolecularBasis& basis) const
{
    CKineticEnergyMatrix kinmat;
    
    if (_locRank == mpi::master())
    {
        // set up GTOs container
        
        CGtoContainer bracontr(molecule, basis);
        
        // compute kinetic energy integrals
        
        kinmat = _compKineticEnergyIntegrals(&bracontr, &bracontr);
    }
    
    return kinmat;
}

CKineticEnergyMatrix
CKineticEnergyIntegralsDriver::compute(const CMolecule&       molecule,
                                       const CMolecularBasis& braBasis,
                                       const CMolecularBasis& ketBasis) const
{
    CKineticEnergyMatrix kinmat;
    
    if (_locRank == mpi::master())
    {
        // set up GTOs containers
        
        CGtoContainer bracontr(molecule, braBasis);
        
        CGtoContainer ketcontr(molecule, ketBasis);
        
        // compute kinetic energy integrals
        
        kinmat = _compKineticEnergyIntegrals(&bracontr, &ketcontr);
    }
    
    return kinmat;
}

CKineticEnergyMatrix
CKineticEnergyIntegralsDriver::compute(const CMolecule&       braMolecule,
                                       const CMolecule&       ketMolecule,
                                       const CMolecularBasis& basis) const
{
    CKineticEnergyMatrix kinmat;
    
    if (_locRank == mpi::master())
    {
        // set up GTOs containers
        
        CGtoContainer bracontr(braMolecule, basis);
        
        CGtoContainer ketcontr(ketMolecule, basis);
        
        // compute kinetic energy integrals
        
        kinmat = _compKineticEnergyIntegrals(&bracontr, &ketcontr);
    }
    
    return kinmat;
}

CKineticEnergyMatrix
CKineticEnergyIntegralsDriver::compute(const CMolecule&       braMolecule,
                                       const CMolecule&       ketMolecule,
                                       const CMolecularBasis& braBasis,
                                       const CMolecularBasis& ketBasis) const
{
    CKineticEnergyMatrix kinmat;
    
    if (_locRank == mpi::master())
    {
        // set up GTOs containers
        
        CGtoContainer bracontr(braMolecule, braBasis);
        
        CGtoContainer ketcontr(ketMolecule, ketBasis);
        
        // compute kinetic energy integrals
        
        kinmat = _compKineticEnergyIntegrals(&bracontr, &ketcontr);
    }
    
    return kinmat;
}

void
CKineticEnergyIntegralsDriver::compute(      double*    intsBatch,
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
    
    // compute kinetic energy integrals
    
    _compKineticEnergyForGtoBlocks(&dist, braGtoBlock, ketGtoBlock); 
}

CKineticEnergyMatrix
CKineticEnergyIntegralsDriver::_compKineticEnergyIntegrals(const CGtoContainer* braGtoContainer,
                                                           const CGtoContainer* ketGtoContainer) const
{
    // check if GTOs containers are same on bra and ket sides
    
    auto symbk = ((*braGtoContainer) == (*ketGtoContainer));
    
    // determine dimensions of overlap matrix
    
    auto nrow = braGtoContainer->getNumberOfAtomicOrbitals();
    
    auto ncol = ketGtoContainer->getNumberOfAtomicOrbitals();
    
    // allocate dense matrix for kinetic energy integrals
    
    CDenseMatrix kinmat(nrow, ncol);
    
    // set up distributio pattern
    
    dist1e dstyp = (symbk) ? dist1e::symsq : dist1e::rect;
    
    COneIntsDistribution* distpat = new COneIntsDistribution(kinmat.values(),
                                                             nrow, ncol, dstyp);
    
    // compute kinetic energy integral blocks
    
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
                        
                        _compKineticEnergyForGtoBlocks(distpat, bgtos, kgtos);
                    }
                }
            }
        }
    }
    
    // deallocate distribution pattern
    
    delete distpat;
    
    return CKineticEnergyMatrix(kinmat);
}

void
CKineticEnergyIntegralsDriver::_compKineticEnergyForGtoBlocks(      COneIntsDistribution* distPattern,
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
    
    CMemBlock2D<double> rfacts(pdim, 4 * pmax);
    
    // set up tensors of PA and PB distances
    
    auto btcomps = intsfunc::getNumberOfComponentsInDistancesTensor(bang);
    
    auto ktcomps = intsfunc::getNumberOfComponentsInDistancesTensor(kang);
    
    auto rpa = (btcomps > 0) ? CMemBlock2D<double>(pdim, btcomps * pmax) : CMemBlock2D<double>();
    
    auto rpb = (ktcomps > 0) ? CMemBlock2D<double>(pdim, ktcomps * pmax) : CMemBlock2D<double>();
    
    // allocate primitives and auxilary integrals buffer
    
    auto ncart = angmom::to_CartesianComponents(bang, kang);
    
    CMemBlock2D<double> primbuffer(pdim, ncart * pmax);
    
    CMemBlock2D<double> auxbuffer(pdim, 2 * pmax);
    
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
        
        intsfunc::compFactorsForKineticEnergy(rfacts, bragtos, ketgtos, i);
        
        // compute tensors of distances: R(PA) = P - A
        
        intsfunc::compTensorsPA(rpa, rab, rfacts, 4, bragtos, ketgtos, i);
        
        // compute tensors of distances: R(PB) = P - B
        
        intsfunc::compTensorsPB(rpb, rab, rfacts, 4, bragtos, ketgtos, i);
        
        // compite primitive kinetic energy integrals
        
        _compPrimKineticEnergyInts(primbuffer, auxbuffer, rfacts, rab, rpa,
                                   rpb, bragtos, ketgtos, i);
        
        // contract primitive kinetic energy integrals
        
        genfunc::contract(cartbuffer, primbuffer, 0, bragtos, ketgtos, i);
        
        // transform Cartesian to spherical integrals
        
        genfunc::transform(spherbuffer, cartbuffer, bmom, kmom, 0, 0, kdim);
        
        // add batch of integrals to integrals matrix
        
        distpat.distribute(spherbuffer, bragtos, ketgtos, symbk, i);
    }
}

void
CKineticEnergyIntegralsDriver::_compPrimKineticEnergyInts(      CMemBlock2D<double>&  primBuffer,
                                                                CMemBlock2D<double>&  auxBuffer,
                                                          const CMemBlock2D<double>&  osFactors,
                                                          const CMemBlock2D<double>&  abDistances,
                                                          const CMemBlock2D<double>&  paDistances,
                                                          const CMemBlock2D<double>&  pbDistances,
                                                          const CGtoBlock&            braGtoBlock,
                                                          const CGtoBlock&            ketGtoBlock,
                                                          const int32_t               iContrGto) const
{
    // set up angular momentum on bra and ket sides
    
    auto bang = braGtoBlock.getAngularMomentum();
    
    auto kang = ketGtoBlock.getAngularMomentum();
    
    // compute (s|T|s) auxilary integrals
    
    kinrecfunc::compKineticEnergyForSS(primBuffer, auxBuffer, osFactors,
                                       abDistances, braGtoBlock, ketGtoBlock,
                                       iContrGto);
    
    // compute (s|T|p) auxilary integrals
    
    if ((bang == 0) && (kang == 1))
    {
       kinrecfunc::compKineticEnergyForSP(primBuffer, auxBuffer, osFactors,
                                          pbDistances, braGtoBlock, ketGtoBlock,
                                          iContrGto);
        
        return;
    }
    
    // compute (s|T|d) auxilary integrals
    
    if ((bang == 0) && (kang == 2))
    {
       kinrecfunc::compKineticEnergyForSD(primBuffer, auxBuffer, osFactors,
                                          pbDistances, braGtoBlock, ketGtoBlock,
                                          iContrGto);
        
        return;
    }
    
    // compute (s|T|f) auxilary integrals
    
    if ((bang == 0) && (kang == 3))
    {
       kinrecfunc::compKineticEnergyForSF(primBuffer, auxBuffer, osFactors,
                                          pbDistances, braGtoBlock, ketGtoBlock,
                                          iContrGto);
        
        return;
    }
    
    // compute (s|T|g) auxilary integrals
    
    if ((bang == 0) && (kang == 4))
    {
       kinrecfunc::compKineticEnergyForSG(primBuffer, auxBuffer, osFactors,
                                          pbDistances, braGtoBlock, ketGtoBlock,
                                          iContrGto);
        
        return;
    }
    
    // compute (p|T|s) auxilary integrals
    
    if ((bang == 1) && (kang == 0))
    {
       kinrecfunc::compKineticEnergyForPS(primBuffer, auxBuffer, osFactors,
                                          paDistances, braGtoBlock, ketGtoBlock,
                                          iContrGto);
        
        return;
    }
    
    // compute (d|T|s) auxilary integrals
    
    if ((bang == 2) && (kang == 0))
    {
       kinrecfunc::compKineticEnergyForDS(primBuffer, auxBuffer, osFactors,
                                          paDistances, braGtoBlock, ketGtoBlock,
                                          iContrGto);
        
        return;
    }
    
    // compute (f|T|s) auxilary integrals
    
    if ((bang == 3) && (kang == 0))
    {
       kinrecfunc::compKineticEnergyForFS(primBuffer, auxBuffer, osFactors,
                                          paDistances, braGtoBlock, ketGtoBlock,
                                          iContrGto);
        
        return;
    }
    
    // compute (g|T|s) auxilary integrals
    
    if ((bang == 4) && (kang == 0))
    {
       kinrecfunc::compKineticEnergyForGS(primBuffer, auxBuffer, osFactors,
                                          paDistances, braGtoBlock, ketGtoBlock,
                                          iContrGto);
        
        return;
    }
    
    // compute (p|T|p) auxilary integrals
    
    if ((bang == 1) && (kang == 1))
    {
        kinrecfunc::compKineticEnergyForPP(primBuffer, auxBuffer, osFactors,
                                           paDistances, pbDistances, braGtoBlock,
                                           ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (p|T|d) auxilary integrals
    
    if ((bang == 1) && (kang == 2))
    {
        kinrecfunc::compKineticEnergyForPD(primBuffer, auxBuffer, osFactors,
                                           paDistances, pbDistances, braGtoBlock,
                                           ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (d|T|p) auxilary integrals
    
    if ((bang == 2) && (kang == 1))
    {
        kinrecfunc::compKineticEnergyForDP(primBuffer, auxBuffer, osFactors,
                                           paDistances, pbDistances, braGtoBlock,
                                           ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (p|T|f) auxilary integrals
    
    if ((bang == 1) && (kang == 3))
    {
        kinrecfunc::compKineticEnergyForPF(primBuffer, auxBuffer, osFactors,
                                           paDistances, pbDistances, braGtoBlock,
                                           ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (f|T|p) auxilary integrals
    
    if ((bang == 3) && (kang == 1))
    {
        kinrecfunc::compKineticEnergyForFP(primBuffer, auxBuffer, osFactors,
                                           paDistances, pbDistances, braGtoBlock,
                                           ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (p|T|g) auxilary integrals
    
    if ((bang == 1) && (kang == 4))
    {
        kinrecfunc::compKineticEnergyForPG(primBuffer, auxBuffer, osFactors,
                                           paDistances, pbDistances, braGtoBlock,
                                           ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (g|T|p) auxilary integrals
    
    if ((bang == 4) && (kang == 1))
    {
        kinrecfunc::compKineticEnergyForGP(primBuffer, auxBuffer, osFactors,
                                           paDistances, pbDistances, braGtoBlock,
                                           ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (d|T|d) auxilary integrals
    
    if ((bang == 2) && (kang == 2))
    {
        kinrecfunc::compKineticEnergyForDD(primBuffer, auxBuffer, osFactors,
                                           paDistances, pbDistances, braGtoBlock,
                                           ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (d|T|f) auxilary integrals
    
    if ((bang == 2) && (kang == 3))
    {
        kinrecfunc::compKineticEnergyForDF(primBuffer, auxBuffer, osFactors,
                                           paDistances, pbDistances, braGtoBlock,
                                           ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (f|T|d) auxilary integrals
    
    if ((bang == 3) && (kang == 2))
    {
        kinrecfunc::compKineticEnergyForFD(primBuffer, auxBuffer, osFactors,
                                           paDistances, pbDistances, braGtoBlock,
                                           ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (d|T|g) auxilary integrals
    
    if ((bang == 2) && (kang == 4))
    {
        kinrecfunc::compKineticEnergyForDG(primBuffer, auxBuffer, osFactors,
                                           paDistances, pbDistances, braGtoBlock,
                                           ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (g|T|d) auxilary integrals
    
    if ((bang == 4) && (kang == 2))
    {
        kinrecfunc::compKineticEnergyForGD(primBuffer, auxBuffer, osFactors,
                                           paDistances, pbDistances, braGtoBlock,
                                           ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (f|T|f) auxilary integrals
    
    if ((bang == 3) && (kang == 3))
    {
        kinrecfunc::compKineticEnergyForFF(primBuffer, auxBuffer, osFactors,
                                           paDistances, pbDistances, braGtoBlock,
                                           ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (f|T|g) auxilary integrals
    
    if ((bang == 3) && (kang == 4))
    {
        kinrecfunc::compKineticEnergyForFG(primBuffer, auxBuffer, osFactors,
                                           paDistances, pbDistances, braGtoBlock,
                                           ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (g|T|f) auxilary integrals
    
    if ((bang == 4) && (kang == 3))
    {
        kinrecfunc::compKineticEnergyForGF(primBuffer, auxBuffer, osFactors,
                                           paDistances, pbDistances, braGtoBlock,
                                           ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (g|T|g) auxilary integrals
    
    if ((bang == 4) && (kang == 4))
    {
        kinrecfunc::compKineticEnergyForGG(primBuffer, auxBuffer, osFactors,
                                           paDistances, pbDistances, braGtoBlock,
                                           ketGtoBlock, iContrGto);
        
        return;
    }
    
    // NOTE: Add l > 4 terms if needed
}

