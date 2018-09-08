//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "ElectronRepulsionIntegralsDriver.hpp"

#include "SystemClock.hpp"
#include "GtoPairsContainer.hpp"
#include "SphericalMomentum.hpp"
#include "GenFunc.hpp"
#include "TwoIntsFunc.hpp"
#include "AngularMomentum.hpp"
#include "EriFuncForSG.hpp"
#include "EriFuncForH.hpp"
#include "EriFuncForI.hpp"
#include "EriFuncForK.hpp"
#include "EriFuncForL.hpp"

CElectronRepulsionIntegralsDriver::CElectronRepulsionIntegralsDriver(const int32_t  globRank,
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

CElectronRepulsionIntegralsDriver::~CElectronRepulsionIntegralsDriver()
{
    
}

void
CElectronRepulsionIntegralsDriver::compute(const CMolecule&       molecule,
                                           const CMolecularBasis& aoBasis,
                                           const double           threshold,
                                                 COutputStream&   oStream,
                                                 MPI_Comm         comm) const
{
    CSystemClock eritim;
    
    // generate GTOs pairs blocks for AO basis on bra side
    
    CGtoPairsContainer bgtopairs(molecule, aoBasis, 1.0e-13);
    
    // split GTOs pairs into batches on bra side
    
    auto bbpairs = bgtopairs.split(5000);
    
    // FIX ME: ....
}

void
CElectronRepulsionIntegralsDriver::compElectronRepulsionForGtoPairsBlocks(const CGtoPairsBlock& braGtoPairsBlock,
                                                                          const CGtoPairsBlock& ketGtoPairsBlock) const
{
    // copy GTOs pairs blocks for bra and ket sides
    
    auto brapairs = braGtoPairsBlock;
    
    auto ketpairs = ketGtoPairsBlock;
    
    // determine symmetry of bra and ket sides
    
    bool symbk = (brapairs == ketpairs);
    
    // set up spherical angular momentum for bra and ket sides
    
    CSphericalMomentum amom(brapairs.getBraAngularMomentum());
    
    CSphericalMomentum bmom(brapairs.getKetAngularMomentum());
    
    CSphericalMomentum cmom(ketpairs.getBraAngularMomentum());
    
    CSphericalMomentum dmom(ketpairs.getKetAngularMomentum());
    
    // allocate prefactors used in Obara-Saika recursion
    
    auto pdim = ketpairs.getNumberOfScreenedPrimPairs();
    
    auto pmax = brapairs.getMaxContractionDepth();
    
    CMemBlock2D<double> rpq(pdim, 3 * pmax);
    
    CMemBlock2D<double> rfacts(pdim, 4 * pmax);
    
    CMemBlock2D<double> rw(pdim, 3 * pmax);
    
    CMemBlock2D<double> rwp(pdim, 3 * pmax);
    
    CMemBlock2D<double> rwq(pdim, 3 * pmax);
    
    // generate horizontal recursion pattern for bra
    
    auto bhrrvec = _getBraHorizontalRecursionPattern(brapairs, ketpairs);
    
    // generate contracted intermediate integrals list
    
    auto tcvec = genfunc::getTriplesFromQuadrupleIndexes(bhrrvec);
    
    // generate horizontal recursion pattern for ket
    
    auto khrrvec = _getKetHorizontalRecursionPattern(tcvec);
    
    // generate primitive intermediate integrals list
    
    auto tpvec = genfunc::getPairsFromTripleIndexes(khrrvec);
    
    // generate vertical recursion pattern
    
    auto vrrvec = _getVerticalRecursionPattern(tpvec);
    
    // set up primitives buffer indexes
    
    std::vector<int32_t> vrridx;
    
    auto nblk = _getIndexesForVerticalRecursionPattern(vrridx, vrrvec, pmax);
    
    // allocate primitives integrals buffer
    
    CMemBlock2D<double> pbuffer(pdim, nblk);
    
    // set up contracted integrals buffer indexes
    
    std::vector<int32_t> tpidx;
    
    nblk = _getIndexesForContractedIntegrals(tpidx, tpvec);
    
    // allocate contracted integrals buffer
    
    auto cdim = ketpairs.getNumberOfScreenedContrPairs();
    
    CMemBlock2D<double> cbuffer(cdim, nblk);
    
    // initialize Boys function evaluator
    
    auto bord = genfunc::maxOrderOfPair(vrrvec, 0, 0);
    
    CBoysFunction bftab(bord);
    
    CMemBlock<double> bargs(pdim);
    
    CMemBlock2D<double> bvals(pdim, bord + 1);
    
    // TESTING INFO:
    
    printf("*** INTEGRALS (%i,%i|%i,%i) Sym %i:\n", amom.getAngularMomentum(),
           bmom.getAngularMomentum(), cmom.getAngularMomentum(),
           dmom.getAngularMomentum(), symbk);

    printf("HRR on bra side:\n");
    
    for (size_t i = 0; i < bhrrvec.size(); i++)
    {
        printf("-> BHRR(%i,%i|%i,%i):\n", bhrrvec[i].first(), bhrrvec[i].second(),
               bhrrvec[i].third(), bhrrvec[i].fourth()); 
    }
    
    printf("Contr. Intermidiates:\n");
    
    for (size_t i = 0; i < tcvec.size(); i++)
    {
        printf("(0,%i|%i,%i):  ", tcvec[i].first(), tcvec[i].second(), tcvec[i].third());
    }
    
    printf("\n");
    
    printf("HRR on ket side:\n");
    
    for (size_t i = 0; i < khrrvec.size(); i++)
    {
        printf("-> KHRR(0,%i|%i,%i):\n", khrrvec[i].first(), khrrvec[i].second(),
               khrrvec[i].third());
    }
    
    printf("Prim. Intermidiates:\n");
    
    for (size_t i = 0; i < tpvec.size(); i++)
    {
        printf("(0,%i|0,%i)^(%i):%i ", tpvec[i].first(), tpvec[i].second(),
               tpvec[i].third(), tpidx[i]);
    }
    
    printf("\n");
    
    printf("VRR:\n");
    
    for (size_t i = 0; i < vrrvec.size(); i++)
    {
        printf("-> VRR(0,%i|0,%i)^(%i): %i pmax %i pdim %i\n", vrrvec[i].first(), vrrvec[i].second(),
               vrrvec[i].third(), vrridx[i], pmax, pdim);
    }
    
    // END OF TESTING INFO
    
    // loop over contracted GTOs ob bra side
    
    for (int32_t i = 0; i < brapairs.getNumberOfScreenedContrPairs(); i++)
    {
        // compute distances: R(PQ) = P - Q
        
        twointsfunc::compDistancesPQ(rpq, brapairs, ketpairs, symbk, i);
        
        // compute Obara-Saika recursion factors
        
        twointsfunc::compFactorsForElectronRepulsion(rfacts, brapairs, ketpairs,
                                                     symbk, i);
        
        // compute coordinates of center W
        
        twointsfunc::compCoordinatesForW(rw, rfacts, 4, brapairs, ketpairs,
                                         symbk, i);
        
        // compute distances: R(WP) = W - P
        
        twointsfunc::compDistancesWP(rwp, rw, brapairs, ketpairs, symbk, i);
        
        // compute distances: R(WQ) = W - Q;
        
        twointsfunc::compDistancesWQ(rwq, rw, brapairs, ketpairs, symbk, i);
        
        // compute primitive electron repulsion integrals
        
        _compPrimElectronRepulsionInts(pbuffer, vrrvec, vrridx, bftab, bargs,
                                       bvals, bord, rfacts, rpq, rwp, rwq,
                                       brapairs, ketpairs, symbk, i);
        
        // contract primitive electron repulsion integrals
        
        genfunc::contract(cbuffer, pbuffer, tpvec, tpidx, vrrvec, vrridx,
                          brapairs, ketpairs, symbk, i);
    }
}

CVecFourIndexes
CElectronRepulsionIntegralsDriver::_getBraHorizontalRecursionPattern(const CGtoPairsBlock& braGtoPairsBlock,
                                                                     const CGtoPairsBlock& ketGtoPairsBlock) const
{
    // set up angular momentum
    
    auto anga = braGtoPairsBlock.getBraAngularMomentum();
    
    auto angb = braGtoPairsBlock.getKetAngularMomentum();
    
    auto angc = ketGtoPairsBlock.getBraAngularMomentum();
    
    auto angd = ketGtoPairsBlock.getKetAngularMomentum();
    
    // set up recursion buffer
    
    CVecFourIndexes recvec;
    
    recvec.reserve((angc + 1) * (angd + 1));
    
    // set up indexing counters
    
    int32_t spos = 0;
    
    int32_t epos = 1;
    
    // set up initial state of recursion buffer
    
    recvec.push_back(CFourIndexes(anga, angb, angc, angd));
    
    while (true)
    {
        // internal new recursion terms counter
        
        int32_t nterms = 0;
        
        // generate bra and ket Obara-Saika recursion terms
        
        for (int32_t i = spos; i < epos; i++)
        {
            CFourIndexes cidx(recvec[i]);
            
            // ((a - 1) b|g(r,r')| c d) term
            
            CFourIndexes t0idx(cidx.first() - 1, cidx.second(),
                                
                               cidx.third(), cidx.fourth());
            
            if (genfunc::addValidAndUniqueQuadruple(recvec, t0idx)) nterms++;
            
            // ((a - 1) (b + 1) |g(r,r')| c d) term
            
            CFourIndexes t1idx(cidx.first() - 1,  cidx.second() + 1,
                                
                               cidx.third(), cidx.fourth());
            
            if (genfunc::addValidAndUniqueQuadruple(recvec, t1idx)) nterms++;
        }
        
        // break loop, all recursion terms are generrated
        
        if (nterms == 0) break;
        
        // update counters
        
        spos  = epos;
        
        epos += nterms;
    }
    
    return recvec;
}

CVecThreeIndexes
CElectronRepulsionIntegralsDriver::_getKetHorizontalRecursionPattern(const CVecThreeIndexes& leadTerms) const
{
    // set up recursion buffer
    
    CVecThreeIndexes recvec = leadTerms;
    
    // set up indexing counters
    
    int32_t spos = 0;
    
    int32_t epos = static_cast<int32_t>(recvec.size());
    
    while (true)
    {
        // internal new recursion terms counter
        
        int32_t nterms = 0;
        
        // generate bra and ket Obara-Saika recursion terms
        
        for (int32_t i = spos; i < epos; i++)
        {
            CThreeIndexes cidx(recvec[i]);
            
            // (0 b |g(r,r')| (c - 1) d) term
            
            CThreeIndexes t0idx(cidx.first(),  cidx.second() - 1,
                                
                                cidx.third());
            
            if (genfunc::addValidAndUniqueTriple(recvec, t0idx)) nterms++;
            
            // (0 b |g(r,r')| (c - 1) (d + 1)) term
            
            CThreeIndexes t1idx(cidx.first(),  cidx.second() - 1,
                                
                                cidx.third() + 1);
            
            if (genfunc::addValidAndUniqueTriple(recvec, t1idx)) nterms++;
        }
        
        // break loop, all recursion terms are generrated
        
        if (nterms == 0) break;
        
        // update counters
        
        spos  = epos;
        
        epos += nterms;
    }
    
    return recvec;
}

CVecThreeIndexes
CElectronRepulsionIntegralsDriver::_getVerticalRecursionPattern(const CVecThreeIndexes& leadTerms) const
{
    // set up recursion buffer
    
    CVecThreeIndexes recvec = leadTerms;
    
    // set up indexing counters
    
    int32_t spos = 0;
    
    int32_t epos = static_cast<int32_t>(recvec.size());
    
    while (true)
    {
        // internal new recursion terms counter
        
        int32_t nterms = 0;
        
        // generate bra and ket Obara-Saika recursion terms
        
        for (int32_t i = spos; i < epos; i++)
        {
            CThreeIndexes cidx(recvec[i]);
            
            if (cidx.first() != 0)
            {
                // (0 (b - 1) |g(r,r')| 0 d)^(m) term
                
                CThreeIndexes t0idx(cidx.first() - 1,  cidx.second(),
                                    
                                    cidx.third());
                
                if (genfunc::addValidAndUniqueTriple(recvec, t0idx)) nterms++;
                
                // (0 (b - 1) |g(r,r')| 0 d)^(m+1) term
                
                CThreeIndexes t1idx(cidx.first() - 1,  cidx.second(),
                                    
                                    cidx.third() + 1);
                
                if (genfunc::addValidAndUniqueTriple(recvec, t1idx)) nterms++;
                
                // (0 (b - 2) |g(r,r')| 0 d)^(m) term
                
                CThreeIndexes t2idx(cidx.first() - 2,  cidx.second(),
                                    
                                    cidx.third());
                
                if (genfunc::addValidAndUniqueTriple(recvec, t2idx)) nterms++;
                
                // (0 (b - 2) |g(r,r')| 0 d)^(m+1) term
                
                CThreeIndexes t3idx(cidx.first() - 2,  cidx.second(),
                                    
                                    cidx.third() + 1);
                
                if (genfunc::addValidAndUniqueTriple(recvec, t3idx)) nterms++;
                
                // (0 (b - 1) |g(r,r')| 0 (d - 1))^(m+1) term
                
                CThreeIndexes tkidx(cidx.first() - 1,  cidx.second() - 1,
                                    
                                    cidx.third() + 1);
                
                if (genfunc::addValidAndUniqueTriple(recvec, tkidx)) nterms++;
            }
            else
            {
                // (0 0 |g(r,r')| 0 (d - 1))^(m) term
                
                CThreeIndexes t0idx(cidx.first(),  cidx.second() - 1,
                                    
                                    cidx.third());
                
                if (genfunc::addValidAndUniqueTriple(recvec, t0idx)) nterms++;
                
                // (0 0 |g(r,r')| 0 (d - 1))^(m+1) term
                
                CThreeIndexes t1idx(cidx.first(),  cidx.second() - 1,
                                    
                                    cidx.third() + 1);
                
                if (genfunc::addValidAndUniqueTriple(recvec, t1idx)) nterms++;
                
                // (0 0 |g(r,r')| 0 (d - 2))^(m) term
                
                CThreeIndexes t2idx(cidx.first(),  cidx.second() - 2,
                                    
                                    cidx.third());
                
                if (genfunc::addValidAndUniqueTriple(recvec, t2idx)) nterms++;
                
                // (0 0 |g(r,r')| 0 (d - 2))^(m+1) term
                
                CThreeIndexes t3idx(cidx.first(),  cidx.second() - 2,
                                    
                                    cidx.third() + 1);
                
                if (genfunc::addValidAndUniqueTriple(recvec, t3idx)) nterms++;
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
CElectronRepulsionIntegralsDriver::_getIndexesForVerticalRecursionPattern(      std::vector<int32_t>& recIndexes,
                                                                          const CVecThreeIndexes&     recPattern,
                                                                          const int32_t               maxPrimPairs) const
{
    // clear vector and reserve memory
    
    recIndexes.clear();
    
    recIndexes.reserve(recPattern.size() + 1);
    
    // loop over recursion pattern
    
    int32_t nblk = 0;
    
    for (size_t i = 0; i < recPattern.size(); i++)
    {
        recIndexes.push_back(nblk);
        
        nblk += maxPrimPairs * angmom::to_CartesianComponents(recPattern[i].first(),
                                                              recPattern[i].second());
    }
    
    return nblk;
}

int32_t
CElectronRepulsionIntegralsDriver::_getIndexesForContractedIntegrals(      std::vector<int32_t>& contrIndexes,
                                                                     const CVecThreeIndexes&     contrListing) const
{
    // clear vector and reserve memory
    
    contrIndexes.clear();
    
    contrIndexes.reserve(contrListing.size() + 1);
    
    // loop over integrals listing
    
    int32_t nblk = 0;
    
    for (size_t i = 0; i < contrListing.size(); i++)
    {
        contrIndexes.push_back(nblk);
        
        nblk += angmom::to_CartesianComponents(contrListing[i].first(),
                                               contrListing[i].second());
    }
    
    return nblk;
}

void
CElectronRepulsionIntegralsDriver::_compPrimElectronRepulsionInts(      CMemBlock2D<double>&  primBuffer,
                                                                  const CVecThreeIndexes&     recPattern,
                                                                  const std::vector<int32_t>& recIndexes,
                                                                  const CBoysFunction&        bfTable,
                                                                        CMemBlock<double>&    bfArguments,
                                                                        CMemBlock2D<double>&  bfValues,
                                                                  const int32_t               bfOrder,
                                                                  const CMemBlock2D<double>&  osFactors,
                                                                  const CMemBlock2D<double>&  pqDistances,
                                                                  const CMemBlock2D<double>&  wpDistances,
                                                                  const CMemBlock2D<double>&  wqDistances,
                                                                  const CGtoPairsBlock&       braGtoPairsBlock,
                                                                  const CGtoPairsBlock&       ketGtoPairsBlock,
                                                                  const bool                  isBraEqualKet,
                                                                  const int32_t               iContrPair) const
{
    // compute (ss|g(r,r')|ss) integrals
    
    erifunc::compElectronRepulsionForSSSS(primBuffer, recPattern, recIndexes,
                                          bfTable, bfArguments, bfValues, bfOrder,
                                          osFactors, pqDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (ss|g(r,r')|sp) integrals
    
    erifunc::compElectronRepulsionForSSSP(primBuffer, recPattern, recIndexes,
                                          wqDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sp|g(r,r')|ss) integrals
    
    erifunc::compElectronRepulsionForSPSS(primBuffer, recPattern, recIndexes,
                                          wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sp|g(r,r')|sp) integrals
    
    erifunc::compElectronRepulsionForSPSP(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (ss|g(r,r')|sd) integrals
    
    erifunc::compElectronRepulsionForSSSD(primBuffer, recPattern, recIndexes,
                                          osFactors, wqDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sd|g(r,r')|ss) integrals
    
    erifunc::compElectronRepulsionForSDSS(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sp|g(r,r')|sd) integrals
    
    erifunc::compElectronRepulsionForSPSD(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sd|g(r,r')|sp) integrals
    
    erifunc::compElectronRepulsionForSDSP(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sd|g(r,r')|sd) integrals
    
    erifunc::compElectronRepulsionForSDSD(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (ss|g(r,r')|sf) integrals
    
    erifunc::compElectronRepulsionForSSSF(primBuffer, recPattern, recIndexes,
                                          osFactors, wqDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sf|g(r,r')|ss) integrals
    
    erifunc::compElectronRepulsionForSFSS(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sp|g(r,r')|sf) integrals
    
    erifunc::compElectronRepulsionForSPSF(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sf|g(r,r')|sp) integrals
    
    erifunc::compElectronRepulsionForSFSP(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sd|g(r,r')|sf) integrals
    
    erifunc::compElectronRepulsionForSDSF(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sf|g(r,r')|sd) integrals
    
    erifunc::compElectronRepulsionForSFSD(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sf|g(r,r')|sf) integrals
    
    erifunc::compElectronRepulsionForSFSF(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (ss|g(r,r')|sg) integrals
    
    erifunc::compElectronRepulsionForSSSG(primBuffer, recPattern, recIndexes,
                                          osFactors, wqDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sg|g(r,r')|ss) integrals
    
    erifunc::compElectronRepulsionForSGSS(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sp|g(r,r')|sg) integrals
    
    erifunc::compElectronRepulsionForSPSG(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sg|g(r,r')|sp) integrals
    
    erifunc::compElectronRepulsionForSGSP(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sd|g(r,r')|sg) integrals
    
    erifunc::compElectronRepulsionForSDSG(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sg|g(r,r')|sd) integrals
    
    erifunc::compElectronRepulsionForSGSD(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sf|g(r,r')|sg) integrals
    
    erifunc::compElectronRepulsionForSFSG(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sg|g(r,r')|sf) integrals
    
    erifunc::compElectronRepulsionForSGSF(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sg|g(r,r')|sg) integrals
    
    erifunc::compElectronRepulsionForSGSG(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (ss|g(r,r')|sh) integrals
    
    erifunc::compElectronRepulsionForSSSH(primBuffer, recPattern, recIndexes,
                                          osFactors, wqDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sh|g(r,r')|ss) integrals
    
    erifunc::compElectronRepulsionForSHSS(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sp|g(r,r')|sh) integrals
    
    erifunc::compElectronRepulsionForSPSH(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sh|g(r,r')|sp) integrals
    
    erifunc::compElectronRepulsionForSHSP(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sd|g(r,r')|sh) integrals
    
    erifunc::compElectronRepulsionForSDSH(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sh|g(r,r')|sd) integrals
    
    erifunc::compElectronRepulsionForSHSD(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sf|g(r,r')|sh) integrals
    
    erifunc::compElectronRepulsionForSFSH(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sh|g(r,r')|sf) integrals
    
    erifunc::compElectronRepulsionForSHSF(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sf|g(r,r')|sh) integrals
    
    erifunc::compElectronRepulsionForSGSH(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sh|g(r,r')|sg) integrals
    
    erifunc::compElectronRepulsionForSHSG(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sh|g(r,r')|sh) integrals
    
    erifunc::compElectronRepulsionForSHSH(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (ss|g(r,r')|si) integrals
    
    erifunc::compElectronRepulsionForSSSI(primBuffer, recPattern, recIndexes,
                                          osFactors, wqDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (si|g(r,r')|ss) integrals
    
    erifunc::compElectronRepulsionForSISS(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sp|g(r,r')|si) integrals
    
    erifunc::compElectronRepulsionForSPSI(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (si|g(r,r')|sp) integrals
    
    erifunc::compElectronRepulsionForSISP(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sd|g(r,r')|si) integrals
    
    erifunc::compElectronRepulsionForSDSI(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (si|g(r,r')|sd) integrals
    
    erifunc::compElectronRepulsionForSISD(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sf|g(r,r')|si) integrals
    
    erifunc::compElectronRepulsionForSFSI(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (si|g(r,r')|sf) integrals
    
    erifunc::compElectronRepulsionForSISF(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sg|g(r,r')|si) integrals
    
    erifunc::compElectronRepulsionForSGSI(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (si|g(r,r')|sg) integrals
    
    erifunc::compElectronRepulsionForSISG(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sh|g(r,r')|si) integrals
    
    erifunc::compElectronRepulsionForSHSI(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (si|g(r,r')|sh) integrals
    
    erifunc::compElectronRepulsionForSISH(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (si|g(r,r')|si) integrals
    
    erifunc::compElectronRepulsionForSISI(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (ss|g(r,r')|sk) integrals
    
    erifunc::compElectronRepulsionForSSSK(primBuffer, recPattern, recIndexes,
                                          osFactors, wqDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sk|g(r,r')|ss) integrals
    
    erifunc::compElectronRepulsionForSKSS(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sp|g(r,r')|sk) integrals
    
    erifunc::compElectronRepulsionForSPSK(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sk|g(r,r')|sp) integrals
    
    erifunc::compElectronRepulsionForSKSP(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sd|g(r,r')|sk) integrals
    
    erifunc::compElectronRepulsionForSDSK(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sk|g(r,r')|sd) integrals
    
    erifunc::compElectronRepulsionForSKSD(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sf|g(r,r')|sk) integrals
    
    erifunc::compElectronRepulsionForSFSK(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sk|g(r,r')|sf) integrals
    
    erifunc::compElectronRepulsionForSKSF(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sg|g(r,r')|sk) integrals
    
    erifunc::compElectronRepulsionForSGSK(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sk|g(r,r')|sg) integrals
    
    erifunc::compElectronRepulsionForSKSG(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sh|g(r,r')|sk) integrals
    
    erifunc::compElectronRepulsionForSHSK(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sk|g(r,r')|sh) integrals
    
    erifunc::compElectronRepulsionForSKSH(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (si|g(r,r')|sk) integrals
    
    erifunc::compElectronRepulsionForSISK(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sk|g(r,r')|si) integrals
    
    erifunc::compElectronRepulsionForSKSI(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sk|g(r,r')|sk) integrals
    
    erifunc::compElectronRepulsionForSKSK(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (ss|g(r,r')|sl) integrals
    
    erifunc::compElectronRepulsionForSSSL(primBuffer, recPattern, recIndexes,
                                          osFactors, wqDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sl|g(r,r')|ss) integrals
    
    erifunc::compElectronRepulsionForSLSS(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sp|g(r,r')|sl) integrals
    
    erifunc::compElectronRepulsionForSPSL(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sl|g(r,r')|sp) integrals
    
    erifunc::compElectronRepulsionForSLSP(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sd|g(r,r')|sl) integrals
    
    erifunc::compElectronRepulsionForSDSL(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sl|g(r,r')|sd) integrals
    
    erifunc::compElectronRepulsionForSLSD(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sf|g(r,r')|sl) integrals
    
    erifunc::compElectronRepulsionForSFSL(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sl|g(r,r')|sf) integrals
    
    erifunc::compElectronRepulsionForSLSF(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sg|g(r,r')|sl) integrals
    
    erifunc::compElectronRepulsionForSGSL(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sl|g(r,r')|sg) integrals
    
    erifunc::compElectronRepulsionForSLSG(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sh|g(r,r')|sl) integrals
    
    erifunc::compElectronRepulsionForSHSL(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sl|g(r,r')|sh) integrals
    
    erifunc::compElectronRepulsionForSLSH(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (si|g(r,r')|sl) integrals
    
    erifunc::compElectronRepulsionForSISL(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sl|g(r,r')|si) integrals
    
    erifunc::compElectronRepulsionForSLSI(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sk|g(r,r')|sl) integrals
    
    erifunc::compElectronRepulsionForSKSL(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sl|g(r,r')|sk) integrals
    
    erifunc::compElectronRepulsionForSLSK(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // compute (sl|g(r,r')|sl) integrals
    
    erifunc::compElectronRepulsionForSLSL(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, isBraEqualKet,
                                          iContrPair);
    
    // add other integrals 
}
