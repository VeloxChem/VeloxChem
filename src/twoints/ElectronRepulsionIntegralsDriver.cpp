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
#include "BoysFunction.hpp"

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
        printf("(0,%i|%i,%i) ", tcvec[i].first(), tcvec[i].second(), tcvec[i].third());
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
        printf("(0,%i|0,%i)^(%i) ", tpvec[i].first(), tpvec[i].second(), tpvec[i].third());
    }
    
    printf("\n");
    
    printf("VRR:\n");
    
    for (size_t i = 0; i < vrrvec.size(); i++)
    {
        printf("-> VRR(0,%i|0,%i)^(%i):\n", vrrvec[i].first(), vrrvec[i].second(),
               vrrvec[i].third());
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
