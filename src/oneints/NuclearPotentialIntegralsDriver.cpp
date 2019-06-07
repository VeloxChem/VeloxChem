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
#include "RecursionFunctionsList.hpp"
#include "TwoCentersRecursionFunctions.hpp"
#include "GenIntsFunc.hpp"

#include "OverlapRecFuncForSX.hpp"
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
    
    // allocate P center coordinates
    
    CMemBlock2D<double> rp(pdim, 3 * pmax);
    
    // set up PA and PB distances
    
    auto rpa = (bang > 0) ? CMemBlock2D<double>(pdim, 3 * pmax) : CMemBlock2D<double>();
    
    auto rpb = (kang > 0) ? CMemBlock2D<double>(pdim, 3 * pmax) : CMemBlock2D<double>();
    
    // set up PC distances
    
    auto rpc = CMemBlock2D<double>(pdim, 3 * pmax);
    
    // allocate primitives and auxilary integrals buffer
    
    auto recmap = _setRecursionMap(bang, kang, pmax);
    
    auto nblock = recmap.getNumberOfComponents();
    
    CMemBlock2D<double> primbuffer(pdim, nblock * pmax);
    
    // allocate primitive integrals accumulation buffer
    
    auto ncart = angmom::to_CartesianComponents(bang, kang);
    
    CMemBlock2D<double> accbuffer(pdim, ncart * pmax);
    
    // allocate contracted Cartesian integrals buffer
    
    auto kdim = ketgtos.getNumberOfContrGtos();
    
    CMemBlock2D<double> cartbuffer(kdim, ncart);
    
    // allocate contracted spherical integrals buffer
    
    auto nspher = angmom::to_SphericalComponents(bang, kang);
    
    CMemBlock2D<double> spherbuffer(kdim, nspher);
    
    // set uo Boys function data
    
    CBoysFunction bftab(bang + kang);
    
    CMemBlock<double> bargs(pdim);
    
    CMemBlock2D<double> bvals(pdim, bang + kang + 1);
    
    // determine bra and ket sides symmetry
    
    bool symbk = (bragtos == ketgtos);
    
    // set up index of primitive integrals
    
    auto pidx = recmap.getIndexOfTerm(gintsfunc::genIntegral({"Nuclear Potential"}, bang, kang, 0));
    
    for (int32_t i = 0; i < bragtos.getNumberOfContrGtos(); i++)
    {
        // compute distances: R(AB) = A - B
        
        intsfunc::compDistancesAB(rab, bragtos, ketgtos, i);
        
        // compute Obara-Saika recursion factors
        
        intsfunc::compFactorsForNuclearPotential(rfacts, bragtos, ketgtos, i);
        
        // compute coordinates of center P
        
        intsfunc::compCoordinatesForP(rp, rfacts, 3, bragtos, ketgtos, i);
        
        // compute distances: R(PA) = P - A
        
        intsfunc::compDistancesPA(rpa, rp, bragtos, ketgtos, i);
        
        // compute distances: R(PB) = P - B
        
        intsfunc::compDistancesPB(rpb, rp, bragtos, ketgtos, i);
        
        // reset primitive integrals accumulation buffer
        
        accbuffer.zero();
        
        // loop over charges
        
        for (int32_t j = 0; j < qvalues.size(); j++)
        {
            // compute distances: R(PC) = P - C
            
            intsfunc::compDistancesPC(rpc, rp, qcoords, bragtos, ketgtos, i, j); 
            
            // compute primitive integrals
            
            _compPrimNuclearPotentialInts(primbuffer, recmap, bftab, bargs, bvals, bang + kang,
                                          rfacts, rab, rpa, rpb, rpc, bragtos, ketgtos, i);
            
            // add scaled contribution to accumulation buffer
            
            _addPointChargeContribution(accbuffer, primbuffer, pidx, qvalues,
                                        bragtos, ketgtos, i, j);
        }
        
        // contract primitive nuclear potential integrals
        
        genfunc::contract(cartbuffer, accbuffer, 0, bragtos, ketgtos, i);
        
        // transform Cartesian to spherical integrals
        
        genfunc::transform(spherbuffer, cartbuffer, bmom, kmom, 0, 0, kdim);
        
        // add batch of integrals to integrals matrix
        
        distpat.distribute(spherbuffer, bragtos, ketgtos, symbk, i);
    }
}

void
CNuclearPotentialIntegralsDriver::_compPrimNuclearPotentialInts(      CMemBlock2D<double>&  primBuffer,
                                                                const CRecursionMap&        recursionMap,
                                                                const CBoysFunction&        bfTable,
                                                                      CMemBlock<double>&    bfArguments,
                                                                      CMemBlock2D<double>&  bfValues,
                                                                const int32_t               bfOrder,
                                                                const CMemBlock2D<double>&  osFactors,
                                                                const CMemBlock2D<double>&  abDistances,
                                                                const CMemBlock2D<double>&  paDistances,
                                                                const CMemBlock2D<double>&  pbDistances,
                                                                const CMemBlock2D<double>&  pcDistances,
                                                                const CGtoBlock&            braGtoBlock,
                                                                const CGtoBlock&            ketGtoBlock,
                                                                const int32_t               iContrGto) const
{
    ovlrecfunc::compOverlapForSS(primBuffer, recursionMap, osFactors, 3, abDistances, braGtoBlock, ketGtoBlock, iContrGto);
    
    npotrecfunc::compNuclearPotentialForSS(primBuffer, recursionMap, bfTable, bfArguments, bfValues, bfOrder,
                                           osFactors, abDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
    
    npotrecfunc::compNuclearPotentialForSP(primBuffer, recursionMap, pbDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
    
    npotrecfunc::compNuclearPotentialForPS(primBuffer, recursionMap, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
    
    npotrecfunc::compNuclearPotentialForSD(primBuffer, recursionMap, osFactors, pbDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
    
    npotrecfunc::compNuclearPotentialForDS(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
    
    npotrecfunc::compNuclearPotentialForSF(primBuffer, recursionMap, osFactors, pbDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
    
    npotrecfunc::compNuclearPotentialForFS(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
    
    npotrecfunc::compNuclearPotentialForSG(primBuffer, recursionMap, osFactors, pbDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
    
    npotrecfunc::compNuclearPotentialForGS(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
    
    npotrecfunc::compNuclearPotentialForPP(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
    
    npotrecfunc::compNuclearPotentialForPD(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
    
    npotrecfunc::compNuclearPotentialForDP(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
    
    npotrecfunc::compNuclearPotentialForPF(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
    
    npotrecfunc::compNuclearPotentialForFP(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
    
    npotrecfunc::compNuclearPotentialForPG(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
    
    npotrecfunc::compNuclearPotentialForGP(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
    
    npotrecfunc::compNuclearPotentialForDD(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
    
    npotrecfunc::compNuclearPotentialForDF(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
    
    npotrecfunc::compNuclearPotentialForFD(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
    
    npotrecfunc::compNuclearPotentialForDG(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
    
    npotrecfunc::compNuclearPotentialForGD(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
    
    npotrecfunc::compNuclearPotentialForFF(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
    
    npotrecfunc::compNuclearPotentialForFG(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
    
    npotrecfunc::compNuclearPotentialForGF(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
    
    npotrecfunc::compNuclearPotentialForGG(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
    
    // NOTE: add l > 4 recursion here
}

void
CNuclearPotentialIntegralsDriver::_addPointChargeContribution(      CMemBlock2D<double>& accBuffer,
                                                              const CMemBlock2D<double>& primBuffer,
                                                              const int32_t              primIndex,
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
            
            auto pbuf = primBuffer.data(primIndex + i * ncart + j);
            
            #pragma omp simd aligned(abuf, pbuf: VLX_ALIGN)
            for (int32_t k = 0; k < nprim; k++)
            {
                abuf[k] += fact * pbuf[k];
            }
        }
    }
}

CRecursionMap
CNuclearPotentialIntegralsDriver::_setRecursionMap(const int32_t braAngularMomentum,
                                                   const int32_t ketAngularMomentum,
                                                   const int32_t maxNumberOfPrimitives) const
{
    CRecursionFunctionsList recfuncs;
    
    recfuncs.add(CRecursionFunction({"Overlap"}, &t2crecfunc::obRecursionForOverlap));
    
    recfuncs.add(CRecursionFunction({"Nuclear Potential"}, &t2crecfunc::obRecursionForNuclearPotential)); 
    
    auto rterm = gintsfunc::genIntegral({"Nuclear Potential"}, braAngularMomentum,
                                        ketAngularMomentum, 0);
    
    return gintsfunc::genRecursionMap(rterm, recblock::cc, maxNumberOfPrimitives,
                                      recfuncs);
}
