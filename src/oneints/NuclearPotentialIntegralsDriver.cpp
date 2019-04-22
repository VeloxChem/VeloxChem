//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "NuclearPotentialIntegralsDriver.hpp"

#include "GenFunc.hpp"
#include "AngularMomentum.hpp"
#include "OneIntsFunc.hpp"
#include "BoysFunction.hpp"
#include "StringFormat.hpp"

#include "NuclearPotentialRecFuncForSX.hpp"
#include "NuclearPotentialRecFuncForPX.hpp"
#include "NuclearPotentialRecFuncForDX.hpp"
#include "NuclearPotentialRecFuncForFF.hpp"
#include "NuclearPotentialRecFuncForFG.hpp"
#include "NuclearPotentialRecFuncForGF.hpp"
#include "NuclearPotentialRecFuncForGG.hpp"

CNuclearPotentialIntegralsDriver::CNuclearPotentialIntegralsDriver(MPI_Comm comm)
{
    _locRank  = mpi::rank(comm);
    
    _locNodes = mpi::nodes(comm);

    mpi::duplicate(comm, &_locComm);
}

CNuclearPotentialIntegralsDriver::~CNuclearPotentialIntegralsDriver()
{
    mpi::destroy(&_locComm);
}

CNuclearPotentialMatrix
CNuclearPotentialIntegralsDriver::compute(const CMolecule&       molecule,
                                          const CMolecularBasis& basis) const
{
    CNuclearPotentialMatrix npotmat;
    
    if (_locRank == mpi::master())
    {
        // set up GTOs container
        
        CGtoContainer bracontr(molecule, basis);
        
        // set up point charges data
        
        auto pcharges = molecule.getCharges();
        
        auto pcoords  = molecule.getCoordinates();
        
        // compute nuclear potential integrals
        
        npotmat = _compNuclearPotentialIntegrals(&pcharges, &pcoords, &bracontr,
                                                 &bracontr);
    }
    
    return npotmat;
}

CNuclearPotentialMatrix
CNuclearPotentialIntegralsDriver::compute(const CMolecule&       molecule,
                                          const CMolecularBasis& basis,
                                          const CMolecule&       pchgMolecule) const
{
    CNuclearPotentialMatrix npotmat;
    
    if (_locRank == mpi::master())
    {
        // set up GTOs container
        
        CGtoContainer bracontr(molecule, basis);
        
        // set up point charges data
        
        auto pcharges = pchgMolecule.getCharges();
        
        auto pcoords  = pchgMolecule.getCoordinates();
        
        // compute nuclear potential integrals
        
        npotmat = _compNuclearPotentialIntegrals(&pcharges, &pcoords, &bracontr,
                                                 &bracontr);
    }
    
    return npotmat;
}

CNuclearPotentialMatrix
CNuclearPotentialIntegralsDriver::compute(const CMolecule&       molecule,
                                          const CMolecularBasis& braBasis,
                                          const CMolecularBasis& ketBasis,
                                          const CMolecule&       pchgMolecule) const
{
    CNuclearPotentialMatrix npotmat;
    
    if (_locRank == mpi::master())
    {
        // set up GTOs container
        
        CGtoContainer bracontr(molecule, braBasis);

        CGtoContainer ketcontr(molecule, ketBasis);
        
        // set up point charges data
        
        auto pcharges = pchgMolecule.getCharges();
        
        auto pcoords  = pchgMolecule.getCoordinates();
        
        // compute nuclear potential integrals
        
        npotmat = _compNuclearPotentialIntegrals(&pcharges, &pcoords, &bracontr,
                                                 &ketcontr);
    }
    
    return npotmat;
}

CNuclearPotentialMatrix
CNuclearPotentialIntegralsDriver::compute(const CMolecule&       braMolecule,
                                          const CMolecule&       ketMolecule,
                                          const CMolecularBasis& basis,
                                          const CMolecule&       pchgMolecule) const
{
    CNuclearPotentialMatrix npotmat;
    
    if (_locRank == mpi::master())
    {
        // set up GTOs container
        
        CGtoContainer bracontr(braMolecule, basis);

        CGtoContainer ketcontr(ketMolecule, basis);
        
        // set up point charges data
        
        auto pcharges = pchgMolecule.getCharges();
        
        auto pcoords  = pchgMolecule.getCoordinates();
        
        // compute nuclear potential integrals
        
        npotmat = _compNuclearPotentialIntegrals(&pcharges, &pcoords, &bracontr,
                                                 &ketcontr);
    }
    
    return npotmat;
}

CNuclearPotentialMatrix
CNuclearPotentialIntegralsDriver::compute(const CMolecule&       braMolecule,
                                          const CMolecule&       ketMolecule,
                                          const CMolecularBasis& braBasis,
                                          const CMolecularBasis& ketBasis,
                                          const CMolecule&       pchgMolecule) const
{
    CNuclearPotentialMatrix npotmat;
    
    if (_locRank == mpi::master())
    {
        // set up GTOs container
        
        CGtoContainer bracontr(braMolecule, braBasis);

        CGtoContainer ketcontr(ketMolecule, ketBasis);
        
        // set up point charges data
        
        auto pcharges = pchgMolecule.getCharges();
        
        auto pcoords  = pchgMolecule.getCoordinates();
        
        // compute nuclear potential integrals
        
        npotmat = _compNuclearPotentialIntegrals(&pcharges, &pcoords, &bracontr,
                                                 &ketcontr);
    }
    
    return npotmat;
}

void
CNuclearPotentialIntegralsDriver::compute(      double*              intsBatch,
                                          const CMemBlock<double>*   charges,
                                          const CMemBlock2D<double>* coordinates,
                                          const CGtoBlock&           braGtoBlock,
                                          const CGtoBlock&           ketGtoBlock) const
{
    // determine dimensions of integrals batch
    
    auto nrow = angmom::to_SphericalComponents(braGtoBlock.getAngularMomentum())
    
              * braGtoBlock.getNumberOfContrGtos();
    
    auto ncol = angmom::to_SphericalComponents(ketGtoBlock.getAngularMomentum())
    
              * ketGtoBlock.getNumberOfContrGtos();
    
    // set up distribution pattern
    
    COneIntsDistribution dist(intsBatch, nrow, ncol, dist1e::batch);
    
    // compute nuclear potential integrals
    
    _compNuclearPotentialForGtoBlocks(&dist, charges, coordinates, braGtoBlock,
                                      ketGtoBlock);
}

CNuclearPotentialMatrix
CNuclearPotentialIntegralsDriver::_compNuclearPotentialIntegrals(const CMemBlock<double>*   charges,
                                                                 const CMemBlock2D<double>* coordinates,
                                                                 const CGtoContainer*       braGtoContainer,
                                                                 const CGtoContainer*       ketGtoContainer) const
{
    // check if GTOs containers are same on bra and ket sides
    
    auto symbk = ((*braGtoContainer) == (*ketGtoContainer));
    
    // determine dimensions of overlap matrix
    
    auto nrow = braGtoContainer->getNumberOfAtomicOrbitals();
    
    auto ncol = ketGtoContainer->getNumberOfAtomicOrbitals();
    
    // allocate dense matrix for nuclear potential integrals
    
    CDenseMatrix npotmat(nrow, ncol);
    
    // set up distributio pattern
    
    dist1e dstyp = (symbk) ? dist1e::symsq : dist1e::rect;
    
    COneIntsDistribution* distpat = new COneIntsDistribution(npotmat.values(),
                                                             nrow, ncol, dstyp);
    
    // compute nuclear potential integral blocks
    
    #pragma omp parallel shared(braGtoContainer, ketGtoContainer, charges,\
                                coordinates, distpat, symbk)
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
                        
                        _compNuclearPotentialForGtoBlocks(distpat, charges, coordinates,
                                                          bgtos, kgtos);
                    }
                }
            }
        }
    }
    
    // deallocate distribution pattern
    
    delete distpat;
    
    return CNuclearPotentialMatrix(npotmat);
}

void
CNuclearPotentialIntegralsDriver::_compNuclearPotentialForGtoBlocks(      COneIntsDistribution* distPattern,
                                                                    const CMemBlock<double>*    charges,
                                                                    const CMemBlock2D<double>*  coordinates,
                                                                    const CGtoBlock&            braGtoBlock,
                                                                    const CGtoBlock&            ketGtoBlock) const
{
    // copy GTOs blocks for bra and ket sides
    
    auto bragtos = braGtoBlock;
    
    auto ketgtos = ketGtoBlock;
    
    // copy distribution pattern
    
    auto distpat = *distPattern;
    
    // copy charges and their coordinates
    
    auto qvalues = *charges;
    
    auto qcoords = *coordinates;
    
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
    
    CMemBlock2D<double> rfacts(pdim, 3 * pmax);
    
    // allocate P coordinates
    
    CMemBlock2D<double> rp(pdim, 3 * pmax);
    
    // set up tensors of PA and PB distances
    
    auto btcomps = intsfunc::getNumberOfComponentsInDistancesTensor(bang);
    
    auto ktcomps = intsfunc::getNumberOfComponentsInDistancesTensor(kang);
    
    auto rpa = (btcomps > 0) ? CMemBlock2D<double>(pdim, btcomps * pmax) : CMemBlock2D<double>();
    
    auto rpb = (ktcomps > 0) ? CMemBlock2D<double>(pdim, ktcomps * pmax) : CMemBlock2D<double>();
    
    auto pcord = kang + bang;
    
    if (pcord == 0) pcord = 1;
    
    auto bkcomps = intsfunc::getNumberOfComponentsInDistancesTensor(pcord);
    
    auto rpc = CMemBlock2D<double>(pdim, bkcomps * pmax);
    
    // Boys function data
    
    auto bord = bang + kang;
    
    CBoysFunction bftab(bord);
    
    CMemBlock<double> bargs(pdim);
    
    CMemBlock2D<double> bvals(pdim, bord + 1);
    
    // allocate primitives and auxilary integrals buffer
    
    auto ncart = angmom::to_CartesianComponents(bang, kang);
    
    CMemBlock2D<double> primbuffer(pdim, ncart * pmax);
    
    CMemBlock2D<double> accbuffer(pdim, ncart * pmax);
    
    CMemBlock2D<double> auxbuffer(pdim, (bord + 1) * pmax);
    
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
        
        intsfunc::compFactorsForNuclearPotential(rfacts, bragtos, ketgtos, i);
        
        // compute coordinates of center P
        
        intsfunc::compCoordinatesForP(rp, rfacts, 3, bragtos, ketgtos, i);
        
        // compute tensors of distances: R(PA) = P - A
        
        intsfunc::compTensorsPA(rpa, rp, bragtos, ketgtos, i);
        
        // compute tnesors of distances: R(PB) = P - B
        
        intsfunc::compTensorsPB(rpb, rp, bragtos, ketgtos, i);
        
        // reset accumulation buffer
        
        accbuffer.zero();
        
        // loop over charges
        
        for (int32_t j = 0; j < qvalues.size(); j++)
        {
            // compute tensors: R(PC) = P - C
            
            intsfunc::compTensorsPC(rpc, rp, qcoords, pcord, bragtos, ketgtos, i, j);
            
            // compute primitive integrals
            
            _compPrimNuclearPotentialInts(primbuffer, auxbuffer, bftab, bargs,
                                          bvals, bord, rfacts, rab, rpa, rpb,
                                          rpc, bkcomps, bragtos, ketgtos, i);
            
            // add scaled contribution to accumulation buffer
            
            _addPointChargeContribution(accbuffer, primbuffer, qvalues,
                                        bragtos, ketgtos, i, j);
        }
        
        // contract primitive overlap integrals
        
        genfunc::contract(cartbuffer, accbuffer, 0, bragtos, ketgtos, i);
        
        // transform Cartesian to spherical integrals
        
        genfunc::transform(spherbuffer, cartbuffer, bmom, kmom, 0, 0, kdim);
        
        // add batch of integrals to integrals matrix
        
        // add batch of integrals to integrals matrix
        
        distpat.distribute(spherbuffer, bragtos, ketgtos, symbk, i);
    }
}

void
CNuclearPotentialIntegralsDriver::_compPrimNuclearPotentialInts(      CMemBlock2D<double>&  primBuffer,
                                                                      CMemBlock2D<double>&  auxBuffer,
                                                                const CBoysFunction&        bfTable,
                                                                      CMemBlock<double>&    bfArguments,
                                                                      CMemBlock2D<double>&  bfValues,
                                                                const int32_t               bfOrder,
                                                                const CMemBlock2D<double>&  osFactors,
                                                                const CMemBlock2D<double>&  abDistances,
                                                                const CMemBlock2D<double>&  paDistances,
                                                                const CMemBlock2D<double>&  pbDistances,
                                                                const CMemBlock2D<double>&  pcDistances,
                                                                const int32_t               pcComponents,
                                                                const CGtoBlock&            braGtoBlock,
                                                                const CGtoBlock&            ketGtoBlock,
                                                                const int32_t               iContrGto) const
{
    // set up angular momentum on bra and ket sides
    
    auto bang = braGtoBlock.getAngularMomentum();
    
    auto kang = ketGtoBlock.getAngularMomentum();
    
    // compute (s|A(0)|s) integrals
    
    npotrecfunc::compNuclearPotentialForSS(primBuffer, auxBuffer, bfTable,
                                           bfArguments, bfValues, bfOrder,
                                           osFactors, abDistances, pcDistances,
                                           pcComponents, braGtoBlock, ketGtoBlock,
                                           iContrGto);
   
    // compute (s|A(0)|p) integrals
   
    if ((bang == 0) && (kang == 1))
    {
        npotrecfunc::compNuclearPotentialForSP(primBuffer, auxBuffer, pbDistances,
                                               pcDistances, braGtoBlock, ketGtoBlock,
                                               iContrGto);
        
        return;
    }
    
    // compute (p|A(0)|s) integrals
    
    if ((bang == 1) && (kang == 0))
    {
        npotrecfunc::compNuclearPotentialForPS(primBuffer, auxBuffer, paDistances,
                                               pcDistances, braGtoBlock, ketGtoBlock,
                                               iContrGto);
        
        return;
    }

    // compute (p|A(0)|p) integrals
    
    if ((bang == 1) && (kang == 1))
    {
        npotrecfunc::compNuclearPotentialForPP(primBuffer, auxBuffer, osFactors,
                                               paDistances, pbDistances, pcDistances,
                                               braGtoBlock, ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (s|A(0)|d) integrals

    if ((bang == 0) && (kang == 2))
    {
        npotrecfunc::compNuclearPotentialForSD(primBuffer, auxBuffer, osFactors,
                                               pbDistances, pcDistances, braGtoBlock,
                                               ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (d|A(0)|s) integrals
   
    if ((bang == 2) && (kang == 0))
    {
        npotrecfunc::compNuclearPotentialForDS(primBuffer, auxBuffer, osFactors,
                                               paDistances, pcDistances, braGtoBlock,
                                               ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (p|A(0)|d) integrals
    
    if ((bang == 1) && (kang == 2))
    {
        npotrecfunc::compNuclearPotentialForPD(primBuffer, auxBuffer, osFactors,
                                               paDistances, pbDistances, pcDistances,
                                               braGtoBlock, ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (d|A(0)|p) integrals
    
    if ((bang == 2) && (kang == 1))
    {
        npotrecfunc::compNuclearPotentialForDP(primBuffer, auxBuffer, osFactors,
                                               paDistances, pbDistances, pcDistances,
                                               braGtoBlock, ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (d|A(0)|d) integrals
    
    if ((bang == 2) && (kang == 2))
    {
        npotrecfunc::compNuclearPotentialForDD(primBuffer, auxBuffer, osFactors,
                                               paDistances, pbDistances, pcDistances,
                                               braGtoBlock, ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (s|A(0)|f) integrals
    
    if ((bang == 0) && (kang == 3))
    {
        npotrecfunc::compNuclearPotentialForSF(primBuffer, auxBuffer, osFactors,
                                               pbDistances, pcDistances, braGtoBlock,
                                               ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (f|A(0)|s) integrals
    
    if ((bang == 3) && (kang == 0))
    {
        npotrecfunc::compNuclearPotentialForFS(primBuffer, auxBuffer, osFactors,
                                               paDistances, pcDistances, braGtoBlock,
                                               ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (p|A(0)|f) integrals
    
    if ((bang == 1) && (kang == 3))
    {
        npotrecfunc::compNuclearPotentialForPF(primBuffer, auxBuffer, osFactors,
                                               paDistances, pbDistances, pcDistances,
                                               braGtoBlock, ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (f|A(0)|p) integrals
    
    if ((bang == 3) && (kang == 1))
    {
        npotrecfunc::compNuclearPotentialForFP(primBuffer, auxBuffer, osFactors,
                                               paDistances, pbDistances, pcDistances,
                                               braGtoBlock, ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (d|A(0)|f) integrals
    
    if ((bang == 2) && (kang == 3))
    {
        npotrecfunc::compNuclearPotentialForDF(primBuffer, auxBuffer, osFactors,
                                               paDistances, pbDistances, pcDistances,
                                               braGtoBlock, ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (f|A(0)|d) integrals
    
    if ((bang == 3) && (kang == 2))
    {
        npotrecfunc::compNuclearPotentialForFD(primBuffer, auxBuffer, osFactors,
                                               paDistances, pbDistances, pcDistances,
                                               braGtoBlock, ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (f|A(0)|f) integrals
    
    if ((bang == 3) && (kang == 3))
    {
        npotrecfunc::compNuclearPotentialForFF(primBuffer, auxBuffer, osFactors,
                                               paDistances, pbDistances, pcDistances,
                                               braGtoBlock, ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (s|A(0)|g) integrals
    
    if ((bang == 0) && (kang == 4))
    {
        npotrecfunc::compNuclearPotentialForSG(primBuffer, auxBuffer, osFactors,
                                               pbDistances, pcDistances, braGtoBlock,
                                               ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (g|A(0)|s) integrals
    
    if ((bang == 4) && (kang == 0))
    {
        npotrecfunc::compNuclearPotentialForGS(primBuffer, auxBuffer, osFactors,
                                               paDistances, pcDistances, braGtoBlock,
                                               ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (p|A(0)|g) integrals
    
    if ((bang == 1) && (kang == 4))
    {
        npotrecfunc::compNuclearPotentialForPG(primBuffer, auxBuffer, osFactors,
                                               paDistances, pbDistances, pcDistances,
                                               braGtoBlock, ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (g|A(0)|p) integrals
    
    if ((bang == 4) && (kang == 1))
    {
        npotrecfunc::compNuclearPotentialForGP(primBuffer, auxBuffer, osFactors,
                                               paDistances, pbDistances, pcDistances,
                                               braGtoBlock, ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (d|A(0)|g) integrals
    
    if ((bang == 2) && (kang == 4))
    {
        npotrecfunc::compNuclearPotentialForDG(primBuffer, auxBuffer, osFactors,
                                               paDistances, pbDistances, pcDistances,
                                               braGtoBlock, ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (g|A(0)|d) integrals
    
    if ((bang == 4) && (kang == 2))
    {
        npotrecfunc::compNuclearPotentialForGD(primBuffer, auxBuffer, osFactors,
                                               paDistances, pbDistances, pcDistances,
                                               braGtoBlock, ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (f|A(0)|g) integrals
    
    if ((bang == 3) && (kang == 4))
    {
        npotrecfunc::compNuclearPotentialForFG(primBuffer, auxBuffer, osFactors,
                                               paDistances, pbDistances, pcDistances,
                                               braGtoBlock, ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (g|A(0)|f) integrals
    
    if ((bang == 4) && (kang == 3))
    {
        npotrecfunc::compNuclearPotentialForGF(primBuffer, auxBuffer, osFactors,
                                               paDistances, pbDistances, pcDistances,
                                               braGtoBlock, ketGtoBlock, iContrGto);
        
        return;
    }
    
    // compute (g|A(0)|g) integrals
    
    if ((bang == 4) && (kang == 4))
    {
        npotrecfunc::compNuclearPotentialForGG(primBuffer, auxBuffer, osFactors,
                                               paDistances, pbDistances, pcDistances,
                                               braGtoBlock, ketGtoBlock, iContrGto);
        
        return;
    }
    
    // NOTE: add l > 4 recursion here
}

void
CNuclearPotentialIntegralsDriver::_addPointChargeContribution(      CMemBlock2D<double>& accBuffer,
                                                              const CMemBlock2D<double>& primBuffer,
                                                              const CMemBlock<double>&   charges,
                                                              const CGtoBlock&           braGtoBlock,
                                                              const CGtoBlock&           ketGtoBlock,
                                                              const int32_t              iContrGto,
                                                              const int32_t              iPointCharge) const
{
    // set up angular momentum for bra and ket sides
    
    auto bang = braGtoBlock.getAngularMomentum();
    
    auto kang = ketGtoBlock.getAngularMomentum();
    
    auto ncart = angmom::to_CartesianComponents(bang, kang);
    
    // set up pointers to primitives data on bra side
    
    auto spos = braGtoBlock.getStartPositions();
    
    auto epos = braGtoBlock.getEndPositions();
    
    auto bdim = epos[iContrGto] - spos[iContrGto];
    
    // set up pointers to primitives data on ket side
    
    auto nprim = ketGtoBlock.getNumberOfPrimGtos();
    
    // set up point charge factor
    
    auto fact = charges.at(iPointCharge);
    
    // loop over dimensions of contracted GTO on bra side
    
    for (int32_t i = 0; i < bdim; i++)
    {
        for (int32_t j = 0; j < ncart; j++)
        {
            auto abuf = accBuffer.data(i * ncart + j);
            
            auto pbuf = primBuffer.data(i * ncart + j);
            
            #pragma omp simd aligned(abuf, pbuf: VLX_ALIGN)
            for (int32_t k = 0; k < nprim; k++)
            {
                abuf[k] += fact * pbuf[k];
            }
        }
    }
}
