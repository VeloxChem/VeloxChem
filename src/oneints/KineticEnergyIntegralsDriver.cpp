//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "KineticEnergyIntegralsDriver.hpp"

#include "GenFunc.hpp"
#include "AngularMomentum.hpp"
#include "OneIntsFunc.hpp"
#include "KineticEnergyRecFunc.hpp"
#include "StringFormat.hpp"

CKineticEnergyIntegralsDriver::CKineticEnergyIntegralsDriver(const int32_t  globRank,
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

CKineticEnergyIntegralsDriver::~CKineticEnergyIntegralsDriver()
{
}

CKineticEnergyMatrix
CKineticEnergyIntegralsDriver::compute(const CMolecule&       molecule,
                                       const CMolecularBasis& basis,
                                             COutputStream&   oStream, 
                                             MPI_Comm         comm) const
{
    CSystemClock timer;
    
    CKineticEnergyMatrix kinmat;
    
    if (_locRank == mpi::master())
    {
        // set up GTOs container
        
        CGtoContainer bracontr(molecule, basis);
        
        // compute kinetic energy integrals
        
        kinmat = _compKineticEnergyIntegrals(&bracontr, &bracontr);
    }
    
    _printComputationTime(timer, oStream);
    
    return kinmat;
}

CKineticEnergyMatrix
CKineticEnergyIntegralsDriver::compute(const CMolecule&       molecule,
                                       const CMolecularBasis& braBasis,
                                       const CMolecularBasis& ketBasis,
                                             COutputStream&   oStream,
                                             MPI_Comm         comm) const
{
    CSystemClock timer;
    
    CKineticEnergyMatrix kinmat;
    
    if (_locRank == mpi::master())
    {
        // set up GTOs containers
        
        CGtoContainer bracontr(molecule, braBasis);
        
        CGtoContainer ketcontr(molecule, ketBasis);
        
        // compute kinetic energy integrals
        
        kinmat = _compKineticEnergyIntegrals(&bracontr, &ketcontr);
    }
    
    _printComputationTime(timer, oStream);
    
    return kinmat;
}

CKineticEnergyMatrix
CKineticEnergyIntegralsDriver::compute(const CMolecule&       braMolecule,
                                       const CMolecule&       ketMolecule,
                                       const CMolecularBasis& basis,
                                             COutputStream&   oStream,
                                             MPI_Comm         comm) const
{
    CSystemClock timer;
    
    CKineticEnergyMatrix kinmat;
    
    if (_locRank == mpi::master())
    {
        // set up GTOs containers
        
        CGtoContainer bracontr(braMolecule, basis);
        
        CGtoContainer ketcontr(ketMolecule, basis);
        
        // compute kinetic energy integrals
        
        kinmat = _compKineticEnergyIntegrals(&bracontr, &ketcontr);
    }
    
    _printComputationTime(timer, oStream);
    
    return kinmat;
}

CKineticEnergyMatrix
CKineticEnergyIntegralsDriver::compute(const CMolecule&       braMolecule,
                                       const CMolecule&       ketMolecule,
                                       const CMolecularBasis& braBasis,
                                       const CMolecularBasis& ketBasis,
                                             COutputStream&   oStream,
                                             MPI_Comm         comm) const
{
    CSystemClock timer;
    
    CKineticEnergyMatrix kinmat;
    
    if (_locRank == mpi::master())
    {
        // set up GTOs containers
        
        CGtoContainer bracontr(braMolecule, braBasis);
        
        CGtoContainer ketcontr(ketMolecule, ketBasis);
        
        // compute kinetic energy integrals
        
        kinmat = _compKineticEnergyIntegrals(&bracontr, &ketcontr);
    }
    
    _printComputationTime(timer, oStream);
    
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
    
    // set up spherical angular momentum for bra and ket sides
    
    CSphericalMomentum bmom(bragtos.getAngularMomentum());
    
    CSphericalMomentum kmom(ketgtos.getAngularMomentum());
    
    // allocate prefactors used in Obara-Saika recursion
    
    auto pdim = ketgtos.getNumberOfPrimGtos();
    
    CMemBlock2D<double> rab(pdim, 3);
    
    auto pmax = bragtos.getMaxContractionDepth();
    
    CMemBlock2D<double> rfacts(pdim, 4 * pmax);
    
    CMemBlock2D<double> rpa(pdim, 3 * pmax);
    
    CMemBlock2D<double> rpb(pdim, 3 * pmax);
    
    // generate recursion pattern
    
    auto recvec = _getRecursionPattern(bragtos, ketgtos);
    
    // set up angular momentum data
    
    auto bang = bragtos.getAngularMomentum();
    
    auto kang = ketgtos.getAngularMomentum();
    
    // set up primitives buffer indexes
    
    std::vector<int32_t> recidx;
    
    auto nblk = _getIndexesForRecursionPattern(recidx, recvec, pmax);
    
    auto pidx = genfunc::findTripleIndex(recidx, recvec, {bang, kang, 0});
        
    // allocate primitives integrals buffer
    
    CMemBlock2D<double> pbuffer(pdim, nblk);
    
    // set up contracted GTOs dimensions
    
    auto kdim = ketgtos.getNumberOfContrGtos();
    
    // allocate contracted Cartesian integrals buffer
    
    auto ncart = angmom::to_CartesianComponents(bang, kang);
    
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
        
        // compute distances: R(PA) = P - A
        
        intsfunc::compDistancesPA(rpa, rab, rfacts, 4, bragtos, ketgtos, i);
        
        // compute distances: R(PB) = P - B
        
        intsfunc::compDistancesPB(rpb, rab, rfacts, 4, bragtos, ketgtos, i);
        
        // compite primitive kinetic energy integrals
        
        _compPrimKineticEnergyInts(pbuffer, recvec, recidx, rfacts, rab, rpa,
                                   rpb, bragtos, ketgtos, i);
        
        // contract primitive kinetic energy integrals
        
        genfunc::contract(cartbuffer, pbuffer, pidx, bragtos, ketgtos, i);
        
        // transform Cartesian to spherical integrals
        
        genfunc::transform(spherbuffer, cartbuffer, bmom, kmom, 0, 0, kdim);
        
        // add batch of integrals to integrals matrix
        
        distpat.distribute(spherbuffer, bragtos, ketgtos, symbk, i);
    }
}

void
CKineticEnergyIntegralsDriver::_compPrimKineticEnergyInts(      CMemBlock2D<double>&  primBuffer,
                                                          const CVecThreeIndexes&     recPattern,
                                                          const std::vector<int32_t>& recIndexes,
                                                          const CMemBlock2D<double>&  osFactors,
                                                          const CMemBlock2D<double>&  abDistances,
                                                          const CMemBlock2D<double>&  paDistances,
                                                          const CMemBlock2D<double>&  pbDistances,
                                                          const CGtoBlock&            braGtoBlock,
                                                          const CGtoBlock&            ketGtoBlock,
                                                          const int32_t               iContrGto) const
{
    // compute (s|t|s) integrals
    
    kinrecfunc::compKineticEnergyForSS(primBuffer, recPattern, recIndexes,
                                       osFactors, abDistances, braGtoBlock,
                                       ketGtoBlock, iContrGto);
    
    // compute (s|p) integrals
    
    kinrecfunc::compKineticEnergyForSP(primBuffer, recPattern, recIndexes,
                                       osFactors, pbDistances, braGtoBlock,
                                       ketGtoBlock, iContrGto);
    
    // compute (p|s) integrals
    
    kinrecfunc::compKineticEnergyForPS(primBuffer, recPattern, recIndexes,
                                       osFactors, paDistances, braGtoBlock,
                                       ketGtoBlock, iContrGto);
    
    // compute (p|p) integrals
    
    kinrecfunc::compKineticEnergyForPP(primBuffer, recPattern, recIndexes,
                                       osFactors, paDistances, braGtoBlock,
                                       ketGtoBlock, iContrGto);
    
    // compute (s|d) integrals
    
    kinrecfunc::compKineticEnergyForSD(primBuffer, recPattern, recIndexes,
                                       osFactors, pbDistances, braGtoBlock,
                                       ketGtoBlock, iContrGto);
    
    // compute (d|s) integrals
    
    kinrecfunc::compKineticEnergyForDS(primBuffer, recPattern, recIndexes,
                                       osFactors, paDistances, braGtoBlock,
                                       ketGtoBlock, iContrGto);
    
    // compute (p|d) integrals
    
    kinrecfunc::compKineticEnergyForPD(primBuffer, recPattern, recIndexes,
                                       osFactors, paDistances, braGtoBlock,
                                       ketGtoBlock, iContrGto);
    
    // compute (d|p) integrals
    
    kinrecfunc::compKineticEnergyForDP(primBuffer, recPattern, recIndexes,
                                       osFactors, paDistances, braGtoBlock,
                                       ketGtoBlock, iContrGto);
    
    // compute (d|d) integrals
    
    kinrecfunc::compKineticEnergyForDD(primBuffer, recPattern, recIndexes,
                                       osFactors, paDistances, braGtoBlock,
                                       ketGtoBlock, iContrGto);
    
    // compute (s|f) integrals
    
    kinrecfunc::compKineticEnergyForSF(primBuffer, recPattern, recIndexes,
                                       osFactors, pbDistances, braGtoBlock,
                                       ketGtoBlock, iContrGto);
    
    // compute (f|s) integrals
    
    kinrecfunc::compKineticEnergyForFS(primBuffer, recPattern, recIndexes,
                                       osFactors, paDistances, braGtoBlock,
                                       ketGtoBlock, iContrGto);
    
    // compute (p|f) integrals
    
    kinrecfunc::compKineticEnergyForPF(primBuffer, recPattern, recIndexes,
                                       osFactors, paDistances, braGtoBlock,
                                       ketGtoBlock, iContrGto);
    
    // compute (f|p) integrals
    
    kinrecfunc::compKineticEnergyForFP(primBuffer, recPattern, recIndexes,
                                       osFactors, paDistances, braGtoBlock,
                                       ketGtoBlock, iContrGto);
    
    // compute (d|f) integrals
    
    kinrecfunc::compKineticEnergyForDF(primBuffer, recPattern, recIndexes,
                                       osFactors, paDistances, braGtoBlock,
                                       ketGtoBlock, iContrGto);
    
    // compute (f|d) integrals
    
    kinrecfunc::compKineticEnergyForFD(primBuffer, recPattern, recIndexes,
                                       osFactors, paDistances, braGtoBlock,
                                       ketGtoBlock, iContrGto);
    
    // compute (f|f) integrals
    
    kinrecfunc::compKineticEnergyForFF(primBuffer, recPattern, recIndexes,
                                       osFactors, paDistances, braGtoBlock,
                                       ketGtoBlock, iContrGto);
    
    // compute (s|g) integrals
    
    kinrecfunc::compKineticEnergyForSG(primBuffer, recPattern, recIndexes,
                                       osFactors, pbDistances, braGtoBlock,
                                       ketGtoBlock, iContrGto);
    
    // compute (g|s) integrals
    
    kinrecfunc::compKineticEnergyForGS(primBuffer, recPattern, recIndexes,
                                       osFactors, paDistances, braGtoBlock,
                                       ketGtoBlock, iContrGto);
    
    // compute (p|g) integrals
    
    kinrecfunc::compKineticEnergyForPG(primBuffer, recPattern, recIndexes,
                                       osFactors, paDistances, braGtoBlock,
                                       ketGtoBlock, iContrGto);
    
    // compute (g|p) integrals
    
    kinrecfunc::compKineticEnergyForGP(primBuffer, recPattern, recIndexes,
                                       osFactors, paDistances, braGtoBlock,
                                       ketGtoBlock, iContrGto);
    
    // compute (d|g) integrals
    
    kinrecfunc::compKineticEnergyForDG(primBuffer, recPattern, recIndexes,
                                       osFactors, paDistances, braGtoBlock,
                                       ketGtoBlock, iContrGto);
    
    // compute (g|d) integrals
    
    kinrecfunc::compKineticEnergyForGD(primBuffer, recPattern, recIndexes,
                                       osFactors, paDistances, braGtoBlock,
                                       ketGtoBlock, iContrGto);
    
    // compute (f|g) integrals
    
    kinrecfunc::compKineticEnergyForFG(primBuffer, recPattern, recIndexes,
                                       osFactors, paDistances, braGtoBlock,
                                       ketGtoBlock, iContrGto);
    
    // compute (g|f) integrals
    
    kinrecfunc::compKineticEnergyForGF(primBuffer, recPattern, recIndexes,
                                       osFactors, paDistances, braGtoBlock,
                                       ketGtoBlock, iContrGto);
    
    // compute (g|g) integrals
    
    kinrecfunc::compKineticEnergyForGG(primBuffer, recPattern, recIndexes,
                                       osFactors, paDistances, braGtoBlock,
                                       ketGtoBlock, iContrGto);
    
    // NOTE: add l > 4 recursion here
}

CVecThreeIndexes
CKineticEnergyIntegralsDriver::_getRecursionPattern(const CGtoBlock& braGtoBlock,
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
            
            if (cidx.third() == 0)
            {
                // kinetic energy recursion
                
                if (cidx.first() != 0)
                {
                    // general recursion for bra and ket sides
                    
                    // (a - 1 |t| b) term
                    
                    CThreeIndexes t1idx(cidx.first() - 1,  cidx.second(), 0);
                    
                    if (genfunc::addValidAndUniqueTriple(recvec, t1idx)) nterms++;
                    
                    // (a - 2 |t| b) term
                    
                    CThreeIndexes t2idx(cidx.first() - 2,  cidx.second(), 0);
                    
                    if (genfunc::addValidAndUniqueTriple(recvec, t2idx)) nterms++;
                    
                    // (a - 1 |t| b - 1) term
                    
                    CThreeIndexes tkidx(cidx.first() - 1,  cidx.second() - 1, 0);
                    
                    if (genfunc::addValidAndUniqueTriple(recvec, tkidx)) nterms++;
                    
                    // (a | b) term
                    
                    CThreeIndexes s0idx(cidx.first(),  cidx.second(), 1);
                    
                    if (genfunc::addValidAndUniqueTriple(recvec, s0idx)) nterms++;
                    
                    // (a - 2 | b) term
                    
                    CThreeIndexes s2idx(cidx.first() - 2,  cidx.second(), 1);
                    
                    if (genfunc::addValidAndUniqueTriple(recvec, s2idx)) nterms++;
                }
                else
                {
                    // reduced recursion for ket side
                    
                    // (0 |t| b - 1) term
                    
                    CThreeIndexes t1idx(cidx.first(),  cidx.second() - 1, 0);
                    
                    if (genfunc::addValidAndUniqueTriple(recvec, t1idx)) nterms++;
                    
                    // (0 |t| b - 2) term
                    
                    CThreeIndexes t2idx(cidx.first(),  cidx.second() - 2, 0);
                    
                    if (genfunc::addValidAndUniqueTriple(recvec, t2idx)) nterms++;
                    
                    // (0 | b) term
                    
                    CThreeIndexes s0idx(cidx.first(),  cidx.second(), 1);
                    
                    if (genfunc::addValidAndUniqueTriple(recvec, s0idx)) nterms++;
                    
                    // (0 | b - 2) term
                    
                    CThreeIndexes s2idx(cidx.first(),  cidx.second() - 2, 1);
                    
                    if (genfunc::addValidAndUniqueTriple(recvec, s2idx)) nterms++;
                }
            }
            else
            {
                // overlap recursion
             
                if (cidx.first() != 0)
                {
                    // general recursion for bra and ket sides
                    
                    // (a - 1 | b) term
                    
                    CThreeIndexes t1idx(cidx.first() - 1,  cidx.second(), 1);
                    
                    if (genfunc::addValidAndUniqueTriple(recvec, t1idx)) nterms++;
                    
                    // (a - 2 | b) term
                    
                    CThreeIndexes t2idx(cidx.first() - 2,  cidx.second(), 1);
                    
                    if (genfunc::addValidAndUniqueTriple(recvec, t2idx)) nterms++;
                    
                    // (a - 1 | b - 1) term
                    
                    CThreeIndexes tkidx(cidx.first() - 1,  cidx.second() - 1, 1);
                    
                    if (genfunc::addValidAndUniqueTriple(recvec, tkidx)) nterms++;
                }
                else
                {
                    // reduced recursion for ket side
                    
                    // (0 | b - 1) term
                    
                    CThreeIndexes t1idx(cidx.first(),  cidx.second() - 1, 1);
                    
                    if (genfunc::addValidAndUniqueTriple(recvec, t1idx)) nterms++;
                    
                    // (0 | b - 2) term
                    
                    CThreeIndexes t2idx(cidx.first(),  cidx.second() - 2, 1);
                    
                    if (genfunc::addValidAndUniqueTriple(recvec, t2idx)) nterms++;
                }
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

int32_t
CKineticEnergyIntegralsDriver::_getIndexesForRecursionPattern(      std::vector<int32_t>& recIndexes,
                                                              const CVecThreeIndexes&     recPattern,
                                                              const int32_t               maxPrimGtos) const
{
    // clear vector and reserve memory
    
    recIndexes.clear();
    
    recIndexes.reserve(recPattern.size() + 1);
    
    // loop over recursion pattern
    
    int32_t nblk = 0;
    
    for (size_t i = 0; i < recPattern.size(); i++)
    {
        recIndexes.push_back(nblk);
        
        nblk += maxPrimGtos * angmom::to_CartesianComponents(recPattern[i].first(),
                                                             recPattern[i].second());
    }
    
    return nblk;
}

void
CKineticEnergyIntegralsDriver::_printComputationTime(const CSystemClock&  timer,
                                                           COutputStream& oStream) const
{
    auto tsec = timer.getElapsedTimeInSeconds();
    
    if (_isLocalMode)
    {
        // FIX ME: we need tags for each driver to be implemented to manage
        //         MPI send/receive cycle.
    }
    
    if (_globRank == mpi::master())
    {
        oStream << fmt::info << "Kinetic energy matrix computed in ";
        
        oStream << fstr::to_string(tsec, 2) << " sec.";
        
        oStream << fmt::end << fmt::blank;
    }
}
