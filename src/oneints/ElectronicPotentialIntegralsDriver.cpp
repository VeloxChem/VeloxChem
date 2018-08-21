//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "ElectronicPotentialIntegralsDriver.hpp"

#include "GenFunc.hpp"
#include "AngularMomentum.hpp"
#include "OneIntsFunc.hpp"
#include "ElectronicPotentialRecFunc.hpp"
#include "StringFormat.hpp"

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
                                                   COutputStream&   oStream, 
                                                   MPI_Comm         comm) const
{
    CSystemClock timer;
    
    CElectronicPotentialMatrix epotmat;
    
    if (_locRank == mpi::master())
    {
        // set up GTOs container
        
        CGtoContainer bracontr(molecule, basis);
        
        // compute electronic potential integrals
        
        epotmat = _compElectronicPotentialIntegrals(&bracontr, &bracontr);
    }
    
    _printComputationTime(timer, oStream);
    
    return epotmat;
}

void
CElectronicPotentialIntegralsDriver::compute(      double*    intsValues,
                                             const CGtoBlock& braGtoBlock,
                                             const CGtoBlock& ketGtoBlock) const
{
    _compElectronicPotentialForGtoBlocks(intsValues, braGtoBlock, ketGtoBlock, 0, 0);
}

CElectronicPotentialMatrix
CElectronicPotentialIntegralsDriver::_compElectronicPotentialIntegrals(const CGtoContainer* braGtoContainer,
                                                                       const CGtoContainer* ketGtoContainer) const
{
    // check if GTOs containers are same on bra and ket sides
    
    auto symbk = ((*braGtoContainer) == (*ketGtoContainer));
    
    // determine dimensions of electronic potential matrix
    
    auto nrow = braGtoContainer->getNumberOfAtomicOrbitals();
    
    auto ncol = ketGtoContainer->getNumberOfAtomicOrbitals();
    
    // allocate dense matrix for electronic potential integrals
    
    CDenseMatrix epotmat(nrow, ncol);
    
    auto epotvals = epotmat.values();
    
    // compute electronic potential integral blocks
    
    #pragma omp parallel shared(braGtoContainer, ketGtoContainer, epotvals,\
                                nrow, ncol, symbk)
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
                        
                        _compElectronicPotentialForGtoBlocks(epotvals, bgtos, kgtos,
                                                             nrow, ncol);
                    }
                }
            }
        }
    }
    
    return CElectronicPotentialMatrix(epotmat);
}

void
CElectronicPotentialIntegralsDriver::_compElectronicPotentialForGtoBlocks(      double*    intsValues,
                                                                          const CGtoBlock& braGtoBlock,
                                                                          const CGtoBlock& ketGtoBlock,
                                                                          const int32_t    nRows,
                                                                          const int32_t    nColumns) const
{
    // copy GTOs blocks for bra and ket sides
    
    auto bragtos = braGtoBlock;
    
    auto ketgtos = ketGtoBlock;
    
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
    
    auto bdim = bragtos.getNumberOfContrGtos();
    
    auto kdim = ketgtos.getNumberOfContrGtos();
    
    // allocate contracted Cartesian integrals buffer
    
    auto ncart = angmom::to_CartesianComponents(bang, kang);
    
    CMemBlock2D<double> cartbuffer(kdim, ncart);
    
    // allocate contracted spherical integrals buffer
    
    auto nspher = angmom::to_SphericalComponents(bang, kang);
    
    CMemBlock2D<double> spherbuffer(kdim, nspher);
    
    // initialize Boys function evaluator
    
    auto bord = genfunc::maxOrderOfPair(recvec, 0, 0);
    
    CBoysFunction bftab(bord);
    
    CMemBlock<double> bargs(pdim);
    
    CMemBlock2D<double> bvals(pdim, bord + 1);
    
    // determine integrals storage scheme
    
    bool diagblk = (bragtos == ketgtos);
    
    bool origord = ((nRows == 0) && (nColumns == 0));
    
    for (int32_t i = 0; i < bragtos.getNumberOfContrGtos(); i++)
    {
        // compute distances: R(AB) = A - B
        
        intsfunc::compDistancesAB(rab, bragtos, ketgtos, i);
        
        // compute Obara-Saika recursion factors
        
        intsfunc::compFactorsForElectronicPotential(rfacts, bragtos, ketgtos, i);
        
        // compute distances: R(PA) = P - A
        
        intsfunc::compDistancesPA(rpa, rab, rfacts, 6, bragtos, ketgtos, i);
        
        // compute distances: R(PB) = P - B
        
        intsfunc::compDistancesPB(rpb, rab, rfacts, 6, bragtos, ketgtos, i);
        
        // compite primitive kinetic energy integrals
        
        _compPrimElectronicPotentialInts(pbuffer, recvec, recidx, bftab, bargs,
                                         bvals, bord, rfacts, rab, rpa, rpb,
                                         bragtos, ketgtos, i);
        
        // contract primitive overlap integrals
        
        genfunc::contract(cartbuffer, pbuffer, pidx, bragtos, ketgtos, i);
        
        // transform Cartesian to spherical integrals
        
        genfunc::transform(spherbuffer, cartbuffer, bmom, kmom, kdim);
        
        // add batch of integrals to integrals matrix
        
        if (origord)
        {
            genfunc::distribute(intsValues, spherbuffer, bang, kang, bdim, kdim, i);
        }
        else
        {
            genfunc::distribute(intsValues, spherbuffer, bragtos, ketgtos,
                                diagblk, nColumns, i);
        }
    }
}

void
CElectronicPotentialIntegralsDriver::_compPrimElectronicPotentialInts(      CMemBlock2D<double>&  primBuffer,
                                                                      const CVecThreeIndexes&     recPattern,
                                                                      const std::vector<int32_t>& recIndexes,
                                                                      const CBoysFunction&        bfTable,
                                                                            CMemBlock<double>&    bfArguments,
                                                                            CMemBlock2D<double>&  bfValues,
                                                                      const int32_t               bfOrder,
                                                                      const CMemBlock2D<double>&  osFactors,
                                                                      const CMemBlock2D<double>&  abDistances,
                                                                      const CMemBlock2D<double>&  paDistances,
                                                                      const CMemBlock2D<double>&  pbDistances,
                                                                      const CGtoBlock&            braGtoBlock,
                                                                      const CGtoBlock&            ketGtoBlock,
                                                                      const int32_t               iContrGto) const
{
    // compute (s|g(r,r')|s) integrals
    
    epotrecfunc::compElectronicPotentialForSS(primBuffer, recPattern, recIndexes,
                                              bfTable, bfArguments, bfValues,
                                              bfOrder, osFactors, abDistances,
                                              braGtoBlock, ketGtoBlock,
                                              iContrGto);
    
    // compute (s|g(r,r')|p) integrals
    
    epotrecfunc::compElectronicPotentialForSP(primBuffer, recPattern, recIndexes,
                                              pbDistances, braGtoBlock, ketGtoBlock,
                                              iContrGto);
    
    // compute (p|g(r,r')|s) integrals
    
    epotrecfunc::compElectronicPotentialForPS(primBuffer, recPattern, recIndexes,
                                              paDistances, braGtoBlock, ketGtoBlock,
                                              iContrGto);
    
    // compute (p|g(r,r')|p) integrals
    
    epotrecfunc::compElectronicPotentialForPP(primBuffer, recPattern, recIndexes,
                                              osFactors, paDistances, braGtoBlock,
                                              ketGtoBlock, iContrGto);
    
    // compute (s|g(r,r')|d) integrals
    
    epotrecfunc::compElectronicPotentialForSD(primBuffer, recPattern, recIndexes,
                                              osFactors, pbDistances, braGtoBlock,
                                              ketGtoBlock, iContrGto);
    
    // compute (d|g(r,r')|s) integrals
    
    epotrecfunc::compElectronicPotentialForDS(primBuffer, recPattern, recIndexes,
                                              osFactors, paDistances, braGtoBlock,
                                              ketGtoBlock, iContrGto);
    
    // compute (p|g(r,r')|d) integrals
    
    epotrecfunc::compElectronicPotentialForPD(primBuffer, recPattern, recIndexes,
                                              osFactors, paDistances, braGtoBlock,
                                              ketGtoBlock, iContrGto);
    
    // compute (d|g(r,r')|p) integrals
    
    epotrecfunc::compElectronicPotentialForDP(primBuffer, recPattern, recIndexes,
                                              osFactors, paDistances, braGtoBlock,
                                              ketGtoBlock, iContrGto);
    
    // compute (d|g(r,r')|d) integrals
    
    epotrecfunc::compElectronicPotentialForDD(primBuffer, recPattern, recIndexes,
                                              osFactors, paDistances, braGtoBlock,
                                              ketGtoBlock, iContrGto);
    
    // compute (s|g(r,r')|f) integrals
    
    epotrecfunc::compElectronicPotentialForSF(primBuffer, recPattern, recIndexes,
                                              osFactors, pbDistances, braGtoBlock,
                                              ketGtoBlock, iContrGto);
    
    // compute (f|g(r,r')|s) integrals
    
    epotrecfunc::compElectronicPotentialForFS(primBuffer, recPattern, recIndexes,
                                              osFactors, paDistances, braGtoBlock,
                                              ketGtoBlock, iContrGto);
    
    // compute (p|g(r,r')|f) integrals
    
    epotrecfunc::compElectronicPotentialForPF(primBuffer, recPattern, recIndexes,
                                              osFactors, paDistances, braGtoBlock,
                                              ketGtoBlock, iContrGto);
    
    // compute (f|g(r,r')|p) integrals
    
    epotrecfunc::compElectronicPotentialForFP(primBuffer, recPattern, recIndexes,
                                              osFactors, paDistances, braGtoBlock,
                                              ketGtoBlock, iContrGto);
    
    // compute (d|g(r,r')|f) integrals
    
    epotrecfunc::compElectronicPotentialForDF(primBuffer, recPattern, recIndexes,
                                              osFactors, paDistances, braGtoBlock,
                                              ketGtoBlock, iContrGto);
    
    // compute (f|g(r,r')|d) integrals
    
    epotrecfunc::compElectronicPotentialForFD(primBuffer, recPattern, recIndexes,
                                              osFactors, paDistances, braGtoBlock,
                                              ketGtoBlock, iContrGto);
    
    // compute (f|g(r,r')|f) integrals
    
    epotrecfunc::compElectronicPotentialForFF(primBuffer, recPattern, recIndexes,
                                              osFactors, paDistances, braGtoBlock,
                                              ketGtoBlock, iContrGto);
    
    // compute (s|g(r,r')|g) integrals
    
    epotrecfunc::compElectronicPotentialForSG(primBuffer, recPattern, recIndexes,
                                              osFactors, pbDistances, braGtoBlock,
                                              ketGtoBlock, iContrGto);
    
    // compute (g|g(r,r')|s) integrals
    
    epotrecfunc::compElectronicPotentialForGS(primBuffer, recPattern, recIndexes,
                                              osFactors, paDistances, braGtoBlock,
                                              ketGtoBlock, iContrGto);
    
    // compute (p|g(r,r')|g) integrals
    
    epotrecfunc::compElectronicPotentialForPG(primBuffer, recPattern, recIndexes,
                                              osFactors, paDistances, braGtoBlock,
                                              ketGtoBlock, iContrGto);
    
    // compute (g|g(r,r')|p) integrals
    
    epotrecfunc::compElectronicPotentialForGP(primBuffer, recPattern, recIndexes,
                                              osFactors, paDistances, braGtoBlock,
                                              ketGtoBlock, iContrGto);
    
    // compute (d|g(r,r')|g) integrals
    
    epotrecfunc::compElectronicPotentialForDG(primBuffer, recPattern, recIndexes,
                                              osFactors, paDistances, braGtoBlock,
                                              ketGtoBlock, iContrGto);
    
    // compute (g|g(r,r')|d) integrals
    
    epotrecfunc::compElectronicPotentialForGD(primBuffer, recPattern, recIndexes,
                                              osFactors, paDistances, braGtoBlock,
                                              ketGtoBlock, iContrGto);
    
    // compute (f|g(r,r')|g) integrals
    
    epotrecfunc::compElectronicPotentialForFG(primBuffer, recPattern, recIndexes,
                                              osFactors, paDistances, braGtoBlock,
                                              ketGtoBlock, iContrGto);
    
    // compute (g|g(r,r')|f) integrals
    
    epotrecfunc::compElectronicPotentialForGF(primBuffer, recPattern, recIndexes,
                                              osFactors, paDistances, braGtoBlock,
                                              ketGtoBlock, iContrGto);
    
    // compute (g|g(r,r')|g) integrals
    
    epotrecfunc::compElectronicPotentialForGG(primBuffer, recPattern, recIndexes,
                                              osFactors, paDistances, braGtoBlock,
                                              ketGtoBlock, iContrGto);
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
                
                // (a - 1 |g(r,r')| b)^(m+1) term
                
                CThreeIndexes t1idx(cidx.first() - 1,  cidx.second(),
                                     
                                    cidx.third() + 1);
                
                if (genfunc::addValidAndUniqueTriple(recvec, t1idx)) nterms++;
                
                // (a - 2 |g(r,r')| b)^(m) term
                
                CThreeIndexes t20idx(cidx.first() - 2,  cidx.second(),
                                     
                                     cidx.third());
                
                if (genfunc::addValidAndUniqueTriple(recvec, t20idx)) nterms++;
                
                // (a - 2 |g(r,r')| b)^(m+1) term
                
                CThreeIndexes t21idx(cidx.first() - 2,  cidx.second(),
                                     
                                     cidx.third() + 1);
                
                if (genfunc::addValidAndUniqueTriple(recvec, t21idx)) nterms++;
                
                // (a - 1 |g(r,r')| b - 1)^(m+1) term
                
                CThreeIndexes tkidx(cidx.first() - 1,  cidx.second() - 1,
                                     
                                     cidx.third() + 1);
                
                if (genfunc::addValidAndUniqueTriple(recvec, tkidx)) nterms++;
            }
            else
            {
                // simplified recursion for ket sides
                
                // (0 |g(r,r')| b - 1)^(m+1) term
                
                CThreeIndexes t1idx(cidx.first(),  cidx.second() - 1,
                                     
                                    cidx.third() + 1);
                
                if (genfunc::addValidAndUniqueTriple(recvec, t1idx)) nterms++;
                
                // (0 |g(r,r')| b - 2)^(m) term
                
                CThreeIndexes t20idx(cidx.first(),  cidx.second() - 2,
                                     
                                     cidx.third());
                
                if (genfunc::addValidAndUniqueTriple(recvec, t20idx)) nterms++;
                
                // (0 |g(r,r')| b - 2)^(m+1) term
                
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

int32_t
CElectronicPotentialIntegralsDriver::_getIndexesForRecursionPattern(      std::vector<int32_t>& recIndexes,
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
CElectronicPotentialIntegralsDriver::_printComputationTime(const CSystemClock&  timer,
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
        oStream << fmt::info << "Electronic potential matrix computed in ";
        
        oStream << fstr::to_string(tsec, 2) << " sec.";
        
        oStream << fmt::end << fmt::blank;
    }
}
