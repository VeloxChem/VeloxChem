//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

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
#include "KetHrrFunc.hpp"
#include "BraHrrFunc.hpp"
#include "StringFormat.hpp"

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
CElectronRepulsionIntegralsDriver::compute(      CAOFockMatrix&       aoFockMatrix,
                                           const CAODensityMatrix&    aoDensityMatrix,
                                           const CMolecule&           molecule,
                                           const CMolecularBasis&     aoBasis,
                                           const CScreeningContainer& screeningContainer,
                                                 COutputStream&       oStream,
                                                 MPI_Comm             comm) const
{
    CSystemClock eritim;
    
    // generate GTOs pairs blocks for AO basis on bra side
    
    CGtoPairsContainer bgtopairs(molecule, aoBasis, 1.0e-13);
    
    // split GTOs pairs into batches on bra side
    
    auto bbpairs = bgtopairs.split();
    
    // compute repulsion integrals
    
    aoFockMatrix.zero(); 
    
    _compElectronRepulsionIntegrals(aoFockMatrix, aoDensityMatrix, &bbpairs,
                                    &bbpairs, &screeningContainer);
    
    aoFockMatrix.symmetrize(); 
}

CScreeningContainer
CElectronRepulsionIntegralsDriver::compute(const ericut           screeningScheme,
                                           const double           threshold,
                                           const CMolecule&       molecule,
                                           const CMolecularBasis& aoBasis,
                                                 COutputStream&   oStream,
                                                 MPI_Comm         comm) const
{
    CSystemClock eritim;
    
    // generate GTOs pairs blocks for AO basis on bra side
    
    CGtoPairsContainer bgtopairs(molecule, aoBasis, 1.0e-13);
    
    // split GTOs pairs into batches on bra side
    
    auto bbpairs = bgtopairs.split();
    
    // allocate temporary buffer for Q values on bra side
    
    auto bqbuff = _getQValuesBuffer(bbpairs);
    
    CVecMemBlock<double> kqbuff;
    
    // compute Q values on bra side
    
    computeMaxQValues(&bqbuff, &kqbuff, &bbpairs, &bbpairs);
    
    // copy Q values from bra to ket side
    
    kqbuff = bqbuff;
    
    // initialize screening container
    
    CScreeningContainer qcont(bqbuff, kqbuff, bbpairs, bbpairs, screeningScheme,
                              threshold);
    
    // print Q values computation timings
    
    _printQValuesTiming(molecule, eritim, oStream); 
    
    return qcont;
}

void
CElectronRepulsionIntegralsDriver::compute(      double*         intsBatch,
                                           const CGtoPairsBlock& braGtoPairsBlock,
                                           const CGtoPairsBlock& ketGtoPairsBlock) const
{
    // set up dimensions of integrals batch
    
    auto nrow = braGtoPairsBlock.getNumberOfScreenedContrPairs();
    
    auto ncol = ketGtoPairsBlock.getNumberOfScreenedContrPairs();
    
    // initialize two electron distributor
    
    CTwoIntsDistribution distpat(intsBatch, nrow, ncol, dist2e::batch);
    
    // set up empty screener
    
    CCauchySchwarzScreener qqdat;
    
    // compute batch of two electron integrals
    
    _compElectronRepulsionForGtoPairsBlocks(distpat, qqdat, braGtoPairsBlock,
                                            ketGtoPairsBlock);
}

void
CElectronRepulsionIntegralsDriver::computeMaxQValues(      CVecMemBlock<double>* braQValuesBuffer,
                                                           CVecMemBlock<double>* ketQValuesBuffer,
                                                     const CGtoPairsContainer*   braGtoPairsContainer,
                                                     const CGtoPairsContainer*   ketGtoPairsContainer) const
{
    // determine symmetry of GTOs pairs containers on bra and ket sides
    
    auto symbk = (*braGtoPairsContainer == *ketGtoPairsContainer);
    
    #pragma omp parallel shared(braGtoPairsContainer, ketGtoPairsContainer)
    {
        #pragma omp single nowait
        {
            // Q values for bra side
            
            auto nbra = braGtoPairsContainer->getNumberOfGtoPairsBlocks();
            
            for (int32_t i = 0; i < nbra; i++)
            {
                #pragma omp task firstprivate(i)
                {
                    auto bqprt = (*braQValuesBuffer)[i].data();
                    
                    auto bpairs = braGtoPairsContainer->getGtoPairsBlock(i);
                    
                    _compMaxQValuesForGtoPairsBlock(bqprt, bpairs);
                }
            }
            
            // Q values for ket side if needed
            
            if (!symbk)
            {
                auto nket = ketGtoPairsContainer->getNumberOfGtoPairsBlocks();
                
                for (int32_t i = 0; i < nket; i++)
                {
                    #pragma omp task firstprivate(i)
                    {
                        auto kqprt = (*ketQValuesBuffer)[i].data();
                        
                        auto kpairs = ketGtoPairsContainer->getGtoPairsBlock(i);
                        
                        _compMaxQValuesForGtoPairsBlock(kqprt, kpairs);
                    }
                }
            }
        }
    }
}

void
CElectronRepulsionIntegralsDriver::_compElectronRepulsionForGtoPairsBlocks(      CTwoIntsDistribution&   distPattern,
                                                                           const CCauchySchwarzScreener& intsScreener,
                                                                           const CGtoPairsBlock&         braGtoPairsBlock,
                                                                           const CGtoPairsBlock&         ketGtoPairsBlock) const
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
    
    // set up horizontal recursion buffer for ket side
    
    std::vector<int32_t> khrridx;
    
    nblk = _getIndexesForKetHRRIntegrals(khrridx, khrrvec);
    
    auto cdim = ketpairs.getNumberOfScreenedContrPairs();
    
    CMemBlock2D<double> khrrbuffer(cdim, nblk);
    
    // initialize R(CD) = C - D distance for horizontal recursion
    
    CMemBlock2D<double> rcd(cdim, 3);
    
    // set up horizontal recursion buffer for bra side
    
    std::vector<int32_t> bhrridx;
    
    nblk = _getIndexesForBraHRRIntegrals(bhrridx, bhrrvec);
    
    CMemBlock2D<double> bhrrbuffer(cdim, nblk);
    
    // initialize R(AB) = A - B distance for horizontal recursion
    
    auto rab = brapairs.getDistancesAB();
    
    // allocate spherical integrals buffer
    
    nblk = angmom::to_SphericalComponents(brapairs.getBraAngularMomentum(),
                                          brapairs.getKetAngularMomentum())
         * angmom::to_SphericalComponents(ketpairs.getBraAngularMomentum(),
                                          ketpairs.getKetAngularMomentum());
    
    CMemBlock2D<double> spherbuffer(cdim, nblk);
    
    // initialize Boys function evaluator
    
    auto bord = genfunc::maxOrderOfPair(vrrvec, 0, 0);
    
    CBoysFunction bftab(bord);
    
    CMemBlock<double> bargs(pdim);
    
    CMemBlock2D<double> bvals(pdim, bord + 1);
    
    // set up integrals screening
    
    bool useqq = !intsScreener.isEmpty();
    
    auto qqpairs = ketpairs;
    
    auto ddpairs = ketpairs;
    
    CMemBlock<int32_t> qqvec(cdim);
    
    CMemBlock<int32_t> qqidx(cdim);
    
    CMemBlock<double> distpq(cdim);
    
    CMemBlock<double> qqden(cdim);
    
    // loop over contracted GTOs ob bra side
    
    for (int32_t i = 0; i < brapairs.getNumberOfScreenedContrPairs(); i++)
    {
        // determine GTOs pairs  effective dimensions on ket side
        
        auto nqpdim = (symbk) ? ketpairs.getNumberOfPrimPairs(i) : pdim;
        
        auto nqcdim = (symbk) ? i + 1 : cdim;
        
        // integrals screening: QQ or QQR scheme 
        
        if (useqq)
        {
            // compute effective distances between GTOs pairs on bra and ket sides
            
            if (intsScreener.getScreeningScheme() == ericut::qqr)
            {
                twointsfunc::compEffectiveDistancesPQ(distpq, brapairs, ketpairs,
                                                      symbk, i);
            }
            
            intsScreener.setScreeningVector(qqvec, distpq, symbk, i);
            
            mathfunc::ordering(qqidx.data(), qqvec.data(), nqcdim);
            
            nqcdim = qqpairs.compress(ketpairs, qqvec, nqcdim);
            
            if (nqcdim > 0) nqpdim = qqpairs.getNumberOfPrimPairs(nqcdim - 1);
        }
        
        // all integrals are vanishing in batch, skip computations
        
        if (nqcdim == 0) continue;
        
        // density based screeing of integrals batches
        
        if (useqq)
        {
            // determine max. density element for GTO pairs on ket side
            
            if ((intsScreener.getScreeningScheme() == ericut::qqden) ||
                (intsScreener.getScreeningScheme() == ericut::qqrden))
            {
                distPattern.getMaxDensityElements(qqden, brapairs, qqpairs,
                                                  symbk, nqcdim, i);
                
                intsScreener.setScreeningVector(qqvec, qqidx, qqden, distpq,
                                                nqcdim, i);
                
                nqcdim = ddpairs.compress(qqpairs, qqvec, nqcdim);
                
                if (nqcdim > 0) nqpdim = ddpairs.getNumberOfPrimPairs(nqcdim - 1);
            }
        }
        
        // all integrals are vanishing in batch, skip computations
        
        if (nqcdim == 0) continue;
        
        // density screened scheme
        
        if ((intsScreener.getScreeningScheme() == ericut::qqden) ||
            (intsScreener.getScreeningScheme() == ericut::qqrden))
        {
            // compute distances: R(PQ) = P - Q
            
            twointsfunc::compDistancesPQ(rpq, brapairs, ddpairs, nqpdim, i);
            
            // compute Obara-Saika recursion factors
            
            twointsfunc::compFactorsForElectronRepulsion(rfacts, brapairs,
                                                         ddpairs, nqpdim, i);
            
            // compute coordinates of center W
            
            twointsfunc::compCoordinatesForW(rw, rfacts, 4, brapairs, ddpairs,
                                             nqpdim, i);
            
            // compute distances: R(WP) = W - P
            
            twointsfunc::compDistancesWP(rwp, rw, brapairs, ddpairs, nqpdim, i);
            
            // compute distances: R(WQ) = W - Q;
            
            twointsfunc::compDistancesWQ(rwq, rw, brapairs, ddpairs, nqpdim, i);
            
            // compute primitive electron repulsion integrals
            
            _compPrimElectronRepulsionInts(pbuffer, vrrvec, vrridx, bftab,
                                           bargs, bvals, bord, rfacts, rpq, rwp,
                                           rwq, brapairs, ddpairs, nqpdim, i);
            
            // contract primitive electron repulsion integrals
            
            genfunc::contract(khrrbuffer, pbuffer, khrrvec, khrridx, vrrvec,
                              vrridx, brapairs, ddpairs, nqpdim, nqcdim, i);
            
            // apply horizontal recursion on ket side
            
            ddpairs.getDistancesAB(rcd, nqcdim);
            
            _applyHRRonKet(khrrbuffer, khrrvec, khrridx, rcd, ddpairs, nqcdim,
                           i);
            
            // transform ket side to spherical form
            
            genfunc::transform_ket(bhrrbuffer, khrrbuffer, cmom, dmom, bhrrvec,
                                   bhrridx, khrrvec, khrridx, ddpairs, nqcdim,
                                   i);
            
            // apply horizontal recursion on bra side
            
            _applyHRRonBra(bhrrbuffer, bhrrvec, bhrridx, rab, ddpairs, nqcdim,
                           i);
            
            // transform bra side to spherical form
            
            genfunc::transform_bra(spherbuffer, bhrrbuffer, amom, bmom, bhrrvec,
                                   bhrridx, ddpairs, nqcdim, i);
            
            // distribute integrals: add distribution or Fock formation code
            
            distPattern.distribute(spherbuffer, brapairs, ddpairs, symbk,
                                   nqcdim, i);
        }
        else
        {
            // compute distances: R(PQ) = P - Q
        
            twointsfunc::compDistancesPQ(rpq, brapairs, qqpairs, nqpdim, i);
        
            // compute Obara-Saika recursion factors
        
            twointsfunc::compFactorsForElectronRepulsion(rfacts, brapairs,
                                                         qqpairs, nqpdim, i);
        
            // compute coordinates of center W
        
            twointsfunc::compCoordinatesForW(rw, rfacts, 4, brapairs, qqpairs,
                                             nqpdim, i);
        
            // compute distances: R(WP) = W - P
        
            twointsfunc::compDistancesWP(rwp, rw, brapairs, qqpairs, nqpdim, i);
        
            // compute distances: R(WQ) = W - Q;
        
            twointsfunc::compDistancesWQ(rwq, rw, brapairs, qqpairs, nqpdim, i);
        
            // compute primitive electron repulsion integrals
        
            _compPrimElectronRepulsionInts(pbuffer, vrrvec, vrridx, bftab,
                                           bargs, bvals, bord, rfacts, rpq, rwp,
                                           rwq, brapairs, qqpairs, nqpdim, i);
        
            // contract primitive electron repulsion integrals
        
            genfunc::contract(khrrbuffer, pbuffer, khrrvec, khrridx, vrrvec,
                              vrridx, brapairs, qqpairs, nqpdim, nqcdim, i);
        
            // apply horizontal recursion on ket side
        
            qqpairs.getDistancesAB(rcd, nqcdim);
        
            _applyHRRonKet(khrrbuffer, khrrvec, khrridx, rcd, qqpairs, nqcdim,
                           i);
        
            // transform ket side to spherical form
        
            genfunc::transform_ket(bhrrbuffer, khrrbuffer, cmom, dmom, bhrrvec,
                                   bhrridx, khrrvec, khrridx, qqpairs, nqcdim,
                                   i);
        
            // apply horizontal recursion on bra side
        
            _applyHRRonBra(bhrrbuffer, bhrrvec, bhrridx, rab, qqpairs, nqcdim,
                           i);
        
            // transform bra side to spherical form
        
            genfunc::transform_bra(spherbuffer, bhrrbuffer, amom, bmom, bhrrvec,
                                   bhrridx, qqpairs, nqcdim, i);
        
            // distribute integrals: add distribution or Fock formation code
        
            distPattern.distribute(spherbuffer, brapairs, qqpairs, symbk,
                                   nqcdim, i);
        }
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
CElectronRepulsionIntegralsDriver::_getIndexesForKetHRRIntegrals(      std::vector<int32_t>& intsIndexes,
                                                                 const CVecThreeIndexes&     intsListing) const
{
    // clear vector and reserve memory
    
    intsIndexes.clear();
    
    intsIndexes.reserve(intsListing.size() + 1);
    
    // loop over integrals listing
    
    int32_t nblk = 0;
    
    for (size_t i = 0; i < intsListing.size(); i++)
    {
        intsIndexes.push_back(nblk);
        
        nblk += angmom::to_CartesianComponents(intsListing[i].first())
        
              * angmom::to_CartesianComponents(intsListing[i].second(),
                                               intsListing[i].third());
    }
    
    return nblk;
}

int32_t
CElectronRepulsionIntegralsDriver::_getIndexesForBraHRRIntegrals(      std::vector<int32_t>& intsIndexes,
                                                                 const CVecFourIndexes&      intsListing) const
{
    // clear vector and reserve memory
    
    intsIndexes.clear();
    
    intsIndexes.reserve(intsListing.size() + 1);
    
    // loop over integrals listing
    
    int32_t nblk = 0;
    
    for (size_t i = 0; i < intsListing.size(); i++)
    {
        intsIndexes.push_back(nblk);
        
        nblk += angmom::to_CartesianComponents(intsListing[i].first(),
                                               intsListing[i].second())
        
              *  angmom::to_SphericalComponents(intsListing[i].third(),
                                                intsListing[i].fourth());
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
                                                                  const int32_t               nKetPrimPairs,
                                                                  const int32_t               iContrPair) const
{
    // compute (ss|g(r,r')|ss) integrals
    
    erifunc::compElectronRepulsionForSSSS(primBuffer, recPattern, recIndexes,
                                          bfTable, bfArguments, bfValues, bfOrder,
                                          osFactors, pqDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (ss|g(r,r')|sp) integrals
    
    erifunc::compElectronRepulsionForSSSP(primBuffer, recPattern, recIndexes,
                                          wqDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sp|g(r,r')|ss) integrals
    
    erifunc::compElectronRepulsionForSPSS(primBuffer, recPattern, recIndexes,
                                          wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sp|g(r,r')|sp) integrals
    
    erifunc::compElectronRepulsionForSPSP(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (ss|g(r,r')|sd) integrals
    
    erifunc::compElectronRepulsionForSSSD(primBuffer, recPattern, recIndexes,
                                          osFactors, wqDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sd|g(r,r')|ss) integrals
    
    erifunc::compElectronRepulsionForSDSS(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sp|g(r,r')|sd) integrals
    
    erifunc::compElectronRepulsionForSPSD(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sd|g(r,r')|sp) integrals
    
    erifunc::compElectronRepulsionForSDSP(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sd|g(r,r')|sd) integrals
    
    erifunc::compElectronRepulsionForSDSD(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (ss|g(r,r')|sf) integrals
    
    erifunc::compElectronRepulsionForSSSF(primBuffer, recPattern, recIndexes,
                                          osFactors, wqDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sf|g(r,r')|ss) integrals
    
    erifunc::compElectronRepulsionForSFSS(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sp|g(r,r')|sf) integrals
    
    erifunc::compElectronRepulsionForSPSF(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sf|g(r,r')|sp) integrals
    
    erifunc::compElectronRepulsionForSFSP(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sd|g(r,r')|sf) integrals
    
    erifunc::compElectronRepulsionForSDSF(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sf|g(r,r')|sd) integrals
    
    erifunc::compElectronRepulsionForSFSD(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sf|g(r,r')|sf) integrals
    
    erifunc::compElectronRepulsionForSFSF(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (ss|g(r,r')|sg) integrals
    
    erifunc::compElectronRepulsionForSSSG(primBuffer, recPattern, recIndexes,
                                          osFactors, wqDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sg|g(r,r')|ss) integrals
    
    erifunc::compElectronRepulsionForSGSS(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sp|g(r,r')|sg) integrals
    
    erifunc::compElectronRepulsionForSPSG(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sg|g(r,r')|sp) integrals
    
    erifunc::compElectronRepulsionForSGSP(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sd|g(r,r')|sg) integrals
    
    erifunc::compElectronRepulsionForSDSG(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sg|g(r,r')|sd) integrals
    
    erifunc::compElectronRepulsionForSGSD(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sf|g(r,r')|sg) integrals
    
    erifunc::compElectronRepulsionForSFSG(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sg|g(r,r')|sf) integrals
    
    erifunc::compElectronRepulsionForSGSF(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sg|g(r,r')|sg) integrals
    
    erifunc::compElectronRepulsionForSGSG(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (ss|g(r,r')|sh) integrals
    
    erifunc::compElectronRepulsionForSSSH(primBuffer, recPattern, recIndexes,
                                          osFactors, wqDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sh|g(r,r')|ss) integrals
    
    erifunc::compElectronRepulsionForSHSS(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sp|g(r,r')|sh) integrals
    
    erifunc::compElectronRepulsionForSPSH(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sh|g(r,r')|sp) integrals
    
    erifunc::compElectronRepulsionForSHSP(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sd|g(r,r')|sh) integrals
    
    erifunc::compElectronRepulsionForSDSH(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sh|g(r,r')|sd) integrals
    
    erifunc::compElectronRepulsionForSHSD(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sf|g(r,r')|sh) integrals
    
    erifunc::compElectronRepulsionForSFSH(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sh|g(r,r')|sf) integrals
    
    erifunc::compElectronRepulsionForSHSF(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sf|g(r,r')|sh) integrals
    
    erifunc::compElectronRepulsionForSGSH(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sh|g(r,r')|sg) integrals
    
    erifunc::compElectronRepulsionForSHSG(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sh|g(r,r')|sh) integrals
    
    erifunc::compElectronRepulsionForSHSH(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (ss|g(r,r')|si) integrals
    
    erifunc::compElectronRepulsionForSSSI(primBuffer, recPattern, recIndexes,
                                          osFactors, wqDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (si|g(r,r')|ss) integrals
    
    erifunc::compElectronRepulsionForSISS(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sp|g(r,r')|si) integrals
    
    erifunc::compElectronRepulsionForSPSI(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (si|g(r,r')|sp) integrals
    
    erifunc::compElectronRepulsionForSISP(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sd|g(r,r')|si) integrals
    
    erifunc::compElectronRepulsionForSDSI(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (si|g(r,r')|sd) integrals
    
    erifunc::compElectronRepulsionForSISD(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sf|g(r,r')|si) integrals
    
    erifunc::compElectronRepulsionForSFSI(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (si|g(r,r')|sf) integrals
    
    erifunc::compElectronRepulsionForSISF(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sg|g(r,r')|si) integrals
    
    erifunc::compElectronRepulsionForSGSI(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (si|g(r,r')|sg) integrals
    
    erifunc::compElectronRepulsionForSISG(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sh|g(r,r')|si) integrals
    
    erifunc::compElectronRepulsionForSHSI(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (si|g(r,r')|sh) integrals
    
    erifunc::compElectronRepulsionForSISH(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (si|g(r,r')|si) integrals
    
    erifunc::compElectronRepulsionForSISI(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (ss|g(r,r')|sk) integrals
    
    erifunc::compElectronRepulsionForSSSK(primBuffer, recPattern, recIndexes,
                                          osFactors, wqDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sk|g(r,r')|ss) integrals
    
    erifunc::compElectronRepulsionForSKSS(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sp|g(r,r')|sk) integrals
    
    erifunc::compElectronRepulsionForSPSK(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sk|g(r,r')|sp) integrals
    
    erifunc::compElectronRepulsionForSKSP(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sd|g(r,r')|sk) integrals
    
    erifunc::compElectronRepulsionForSDSK(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sk|g(r,r')|sd) integrals
    
    erifunc::compElectronRepulsionForSKSD(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sf|g(r,r')|sk) integrals
    
    erifunc::compElectronRepulsionForSFSK(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sk|g(r,r')|sf) integrals
    
    erifunc::compElectronRepulsionForSKSF(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sg|g(r,r')|sk) integrals
    
    erifunc::compElectronRepulsionForSGSK(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sk|g(r,r')|sg) integrals
    
    erifunc::compElectronRepulsionForSKSG(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sh|g(r,r')|sk) integrals
    
    erifunc::compElectronRepulsionForSHSK(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sk|g(r,r')|sh) integrals
    
    erifunc::compElectronRepulsionForSKSH(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (si|g(r,r')|sk) integrals
    
    erifunc::compElectronRepulsionForSISK(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sk|g(r,r')|si) integrals
    
    erifunc::compElectronRepulsionForSKSI(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sk|g(r,r')|sk) integrals
    
    erifunc::compElectronRepulsionForSKSK(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (ss|g(r,r')|sl) integrals
    
    erifunc::compElectronRepulsionForSSSL(primBuffer, recPattern, recIndexes,
                                          osFactors, wqDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sl|g(r,r')|ss) integrals
    
    erifunc::compElectronRepulsionForSLSS(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sp|g(r,r')|sl) integrals
    
    erifunc::compElectronRepulsionForSPSL(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sl|g(r,r')|sp) integrals
    
    erifunc::compElectronRepulsionForSLSP(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sd|g(r,r')|sl) integrals
    
    erifunc::compElectronRepulsionForSDSL(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sl|g(r,r')|sd) integrals
    
    erifunc::compElectronRepulsionForSLSD(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sf|g(r,r')|sl) integrals
    
    erifunc::compElectronRepulsionForSFSL(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sl|g(r,r')|sf) integrals
    
    erifunc::compElectronRepulsionForSLSF(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sg|g(r,r')|sl) integrals
    
    erifunc::compElectronRepulsionForSGSL(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sl|g(r,r')|sg) integrals
    
    erifunc::compElectronRepulsionForSLSG(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sh|g(r,r')|sl) integrals
    
    erifunc::compElectronRepulsionForSHSL(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sl|g(r,r')|sh) integrals
    
    erifunc::compElectronRepulsionForSLSH(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (si|g(r,r')|sl) integrals
    
    erifunc::compElectronRepulsionForSISL(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sl|g(r,r')|si) integrals
    
    erifunc::compElectronRepulsionForSLSI(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sk|g(r,r')|sl) integrals
    
    erifunc::compElectronRepulsionForSKSL(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sl|g(r,r')|sk) integrals
    
    erifunc::compElectronRepulsionForSLSK(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // compute (sl|g(r,r')|sl) integrals
    
    erifunc::compElectronRepulsionForSLSL(primBuffer, recPattern, recIndexes,
                                          osFactors, wpDistances, braGtoPairsBlock,
                                          ketGtoPairsBlock, nKetPrimPairs,
                                          iContrPair);
    
    // add other integrals 
}

void
CElectronRepulsionIntegralsDriver::_applyHRRonKet(      CMemBlock2D<double>&  ketBuffer,
                                                  const CVecThreeIndexes&     recPattern,
                                                  const std::vector<int32_t>& recIndexes,
                                                  const CMemBlock2D<double>&  cdDistances,
                                                  const CGtoPairsBlock&       ketGtoPairsBlock,
                                                  const int32_t               nKetContrPairs,
                                                  const int32_t               iContrPair) const
{
    // compute (sx|g(r,r')|pp) integrals
    
    kethrrfunc::compElectronRepulsionForSXPP(ketBuffer, recPattern, recIndexes,
                                             cdDistances, ketGtoPairsBlock,
                                             nKetContrPairs, iContrPair);
    
    // compute (sx|g(r,r')|pd) integrals
    
    kethrrfunc::compElectronRepulsionForSXPD(ketBuffer, recPattern, recIndexes,
                                             cdDistances, ketGtoPairsBlock,
                                             nKetContrPairs, iContrPair);
    
    // compute (sx|g(r,r')|pf) integrals
    
    kethrrfunc::compElectronRepulsionForSXPF(ketBuffer, recPattern, recIndexes,
                                             cdDistances, ketGtoPairsBlock,
                                             nKetContrPairs, iContrPair);
    
    // compute (sx|g(r,r')|pg) integrals
    
    kethrrfunc::compElectronRepulsionForSXPG(ketBuffer, recPattern, recIndexes,
                                             cdDistances, ketGtoPairsBlock,
                                             nKetContrPairs, iContrPair);
    
    // compute (sx|g(r,r')|ph) integrals
    
    kethrrfunc::compElectronRepulsionForSXPH(ketBuffer, recPattern, recIndexes,
                                             cdDistances, ketGtoPairsBlock,
                                             nKetContrPairs, iContrPair);
    
    // compute (sx|g(r,r')|pi) integrals
    
    kethrrfunc::compElectronRepulsionForSXPI(ketBuffer, recPattern, recIndexes,
                                             cdDistances, ketGtoPairsBlock,
                                             nKetContrPairs, iContrPair);
    
    // compute (sx|g(r,r')|pk) integrals
    
    kethrrfunc::compElectronRepulsionForSXPK(ketBuffer, recPattern, recIndexes,
                                             cdDistances, ketGtoPairsBlock,
                                             nKetContrPairs, iContrPair);
    
    // compute (sx|g(r,r')|dd) integrals
    
    kethrrfunc::compElectronRepulsionForSXDD(ketBuffer, recPattern, recIndexes,
                                             cdDistances, ketGtoPairsBlock,
                                             nKetContrPairs, iContrPair);
    
    // compute (sx|g(r,r')|df) integrals
    
    kethrrfunc::compElectronRepulsionForSXDF(ketBuffer, recPattern, recIndexes,
                                             cdDistances, ketGtoPairsBlock,
                                             nKetContrPairs, iContrPair);
    
    // compute (sx|g(r,r')|dg) integrals
    
    kethrrfunc::compElectronRepulsionForSXDG(ketBuffer, recPattern, recIndexes,
                                             cdDistances, ketGtoPairsBlock,
                                             nKetContrPairs, iContrPair);
    
    // compute (sx|g(r,r')|dh) integrals
    
    kethrrfunc::compElectronRepulsionForSXDH(ketBuffer, recPattern, recIndexes,
                                             cdDistances, ketGtoPairsBlock,
                                             nKetContrPairs, iContrPair);
    
    // compute (sx|g(r,r')|di) integrals
    
    kethrrfunc::compElectronRepulsionForSXDI(ketBuffer, recPattern, recIndexes,
                                             cdDistances, ketGtoPairsBlock,
                                             nKetContrPairs, iContrPair);
    
    // compute (sx|g(r,r')|ff) integrals
    
    kethrrfunc::compElectronRepulsionForSXFF(ketBuffer, recPattern, recIndexes,
                                             cdDistances, ketGtoPairsBlock,
                                             nKetContrPairs, iContrPair);
    
    // compute (sx|g(r,r')|fg) integrals
    
    kethrrfunc::compElectronRepulsionForSXFG(ketBuffer, recPattern, recIndexes,
                                             cdDistances, ketGtoPairsBlock,
                                             nKetContrPairs, iContrPair);
    
    // compute (sx|g(r,r')|fh) integrals
    
    kethrrfunc::compElectronRepulsionForSXFH(ketBuffer, recPattern, recIndexes,
                                             cdDistances, ketGtoPairsBlock,
                                             nKetContrPairs, iContrPair);
    
    // compute (sx|g(r,r')|gg) integrals
    
    kethrrfunc::compElectronRepulsionForSXGG(ketBuffer, recPattern, recIndexes,
                                             cdDistances, ketGtoPairsBlock,
                                             nKetContrPairs, iContrPair);
}

void
CElectronRepulsionIntegralsDriver::_applyHRRonBra(      CMemBlock2D<double>&  braBuffer,
                                                  const CVecFourIndexes&      recPattern,
                                                  const std::vector<int32_t>& recIndexes,
                                                  const CMemBlock2D<double>&  abDistances,
                                                  const CGtoPairsBlock&       ketGtoPairsBlock,
                                                  const int32_t               nKetContrPairs,
                                                  const int32_t               iContrPair) const
{
    // compute (pp|g(r,r')|xx) integrals
    
    brahrrfunc::compElectronRepulsionForPPXX(braBuffer, recPattern, recIndexes,
                                             abDistances, ketGtoPairsBlock,
                                             nKetContrPairs, iContrPair);
    
    // compute (pd|g(r,r')|xx) integrals
    
    brahrrfunc::compElectronRepulsionForPDXX(braBuffer, recPattern, recIndexes,
                                             abDistances, ketGtoPairsBlock,
                                             nKetContrPairs, iContrPair);
    
    // compute (pf|g(r,r')|xx) integrals
    
    brahrrfunc::compElectronRepulsionForPFXX(braBuffer, recPattern, recIndexes,
                                             abDistances, ketGtoPairsBlock,
                                             nKetContrPairs, iContrPair);
    
    // compute (pg|g(r,r')|xx) integrals
    
    brahrrfunc::compElectronRepulsionForPGXX(braBuffer, recPattern, recIndexes,
                                             abDistances, ketGtoPairsBlock,
                                             nKetContrPairs, iContrPair);
    
    // compute (ph|g(r,r')|xx) integrals
    
    brahrrfunc::compElectronRepulsionForPHXX(braBuffer, recPattern, recIndexes,
                                             abDistances, ketGtoPairsBlock,
                                             nKetContrPairs, iContrPair);
    
    // compute (pi|g(r,r')|xx) integrals
    
    brahrrfunc::compElectronRepulsionForPIXX(braBuffer, recPattern, recIndexes,
                                             abDistances, ketGtoPairsBlock,
                                             nKetContrPairs, iContrPair);
    
    // compute (pk|g(r,r')|xx) integrals
    
    brahrrfunc::compElectronRepulsionForPKXX(braBuffer, recPattern, recIndexes,
                                             abDistances, ketGtoPairsBlock,
                                             nKetContrPairs, iContrPair);
    
    // compute (dd|g(r,r')|xx) integrals
    
    brahrrfunc::compElectronRepulsionForDDXX(braBuffer, recPattern, recIndexes,
                                             abDistances, ketGtoPairsBlock,
                                             nKetContrPairs, iContrPair);
    
    // compute (df|g(r,r')|xx) integrals
    
    brahrrfunc::compElectronRepulsionForDFXX(braBuffer, recPattern, recIndexes,
                                             abDistances, ketGtoPairsBlock,
                                             nKetContrPairs, iContrPair);
    
    // compute (dg|g(r,r')|xx) integrals
    
    brahrrfunc::compElectronRepulsionForDGXX(braBuffer, recPattern, recIndexes,
                                             abDistances, ketGtoPairsBlock,
                                             nKetContrPairs, iContrPair);
    
    // compute (dh|g(r,r')|xx) integrals
    
    brahrrfunc::compElectronRepulsionForDHXX(braBuffer, recPattern, recIndexes,
                                             abDistances, ketGtoPairsBlock,
                                             nKetContrPairs, iContrPair);
    
    // compute (di|g(r,r')|xx) integrals
    
    brahrrfunc::compElectronRepulsionForDIXX(braBuffer, recPattern, recIndexes,
                                             abDistances, ketGtoPairsBlock,
                                             nKetContrPairs, iContrPair);
    
    // compute (ff|g(r,r')|xx) integrals
    
    brahrrfunc::compElectronRepulsionForFFXX(braBuffer, recPattern, recIndexes,
                                             abDistances, ketGtoPairsBlock,
                                             nKetContrPairs, iContrPair);
    
    // compute (fg|g(r,r')|xx) integrals
    
    brahrrfunc::compElectronRepulsionForFGXX(braBuffer, recPattern, recIndexes,
                                             abDistances, ketGtoPairsBlock,
                                             nKetContrPairs, iContrPair);
    
    // compute (fh|g(r,r')|xx) integrals
    
    brahrrfunc::compElectronRepulsionForFHXX(braBuffer, recPattern, recIndexes,
                                             abDistances, ketGtoPairsBlock,
                                             nKetContrPairs, iContrPair);
    
    // compute (gg|g(r,r')|xx) integrals
    
    brahrrfunc::compElectronRepulsionForGGXX(braBuffer, recPattern, recIndexes,
                                             abDistances, ketGtoPairsBlock,
                                             nKetContrPairs, iContrPair);
}

void
CElectronRepulsionIntegralsDriver::_startHeader(const CGtoPairsContainer& gtoPairs,
                                                      COutputStream&      oStream) const
{
    oStream << fmt::header << "Electron Repulsion Integrals" << fmt::end;
    
    oStream << std::string(27, '=') << fmt::end << fmt::blank;
    
    // GTO pairs screening information
    
    gtoPairs.printScreeningInfo(oStream);
    
    oStream << fmt::blank;
}

void
CElectronRepulsionIntegralsDriver::_printTiming(const CMolecule&     molecule,
                                                const CSystemClock&  timer,
                                                      COutputStream& oStream) const
{
    // NOTE: Silent for local execution mode
    
    if (_isLocalMode) return;
    
    // collect timing data from MPI nodes
    
    auto tsec = timer.getElapsedTimeInSeconds();
    
    CMemBlock<double> tvec;
    
    if (_globRank == mpi::master()) tvec = CMemBlock<double>(_globNodes);
    
    mpi::gather(tvec.data(), tsec, _globRank, MPI_COMM_WORLD);
    
    // print timing data
    
    if (_globRank == mpi::master())
    {
        auto natoms = molecule.getNumberOfAtoms();
        
        std::string str("Two-Electron Integrals Evaluation Timings: ");
        
        oStream << fstr::format(str, 80, fmt::left) << fmt::end << fmt::blank;
        
        for (int32_t i = 0; i < _globNodes; i++)
        {
            // node information
            
            str.assign("MPI Node: ");
            
            str.append(fstr::to_string(i, 3, fmt::left));
            
            // atom batches information
            
            auto nodatm = mpi::batch_size(natoms, i, _globNodes);
            
            auto nodoff = mpi::batch_offset(natoms, i, _globNodes);
            
            str.append(" Atoms in batch: ");
            
            std::string bstr(std::to_string(nodoff));
            
            bstr.append("-");
            
            bstr.append(std::to_string(nodoff + nodatm));
            
            str.append(fstr::format(bstr, 8, fmt::left));
            
            // evaluation time info
            
            str.append(" Time: ");
            
            str.append(fstr::to_string(tvec.at(i), 2));
            
            str.append(" sec.");
            
            oStream << fstr::format(str, 80, fmt::left) << fmt::end;
        }
    }
}

void
CElectronRepulsionIntegralsDriver::_printFockTiming(const CAOFockMatrix& fockMatrix,
                                                    const CSystemClock&  timer,
                                                          COutputStream& oStream) const
{
    // NOTE: Silent for local execution mode
    
    if (_isLocalMode) return;
    
    // collect timing data from MPI nodes
    
    auto tsec = timer.getElapsedTimeInSeconds();
    
    // print timing data
    
    if (_globRank == mpi::master())
    {
        std::string str("AO Fock timing: ");
        
        auto nfock = fockMatrix.getNumberOfFockMatrices();
        
        str.append(std::to_string(nfock));
        
        if (nfock == 1)
        {
            str.append(" matrix ");
        }
        else
        {
            str.append(" matrices ");
        }
        
        str.append("computed in ");
        
        str.append(fstr::to_string(tsec, 2));
        
        str.append(" sec.");
        
        oStream << fstr::format(str, 80, fmt::left) << fmt::end;
    }
}

void
CElectronRepulsionIntegralsDriver::_printQValuesTiming(const CMolecule&     molecule,
                                                       const CSystemClock&  timer,
                                                             COutputStream& oStream) const
{
    // NOTE: Silent for local execution mode
    
    if (_isLocalMode) return;
    
    // collect timing data from MPI nodes
    
    if (_globRank == mpi::master())
    {
        auto tsec = timer.getElapsedTimeInSeconds();
        
        oStream << fmt::info << "ERIs screening factors computed in ";
        
        oStream << fstr::to_string(tsec, 2) << " sec.";
        
        oStream << fmt::end << fmt::blank;
    }
}

void
CElectronRepulsionIntegralsDriver::_compElectronRepulsionIntegrals(      CAOFockMatrix&       aoFockMatrix,
                                                                   const CAODensityMatrix&    aoDensityMatrix,
                                                                   const CGtoPairsContainer*  braGtoPairsContainer,
                                                                   const CGtoPairsContainer*  ketGtoPairsContainer,
                                                                   const CScreeningContainer* screeningContainer) const
{
    // set up pointers to AO Fock and AO density matrices
    
    auto pfock = &aoFockMatrix;
    
    auto pden = &aoDensityMatrix;
    
    #pragma omp parallel shared(braGtoPairsContainer, ketGtoPairsContainer, pfock, pden)
    {
        #pragma omp single nowait
        {
            // determine number of GTOs pairs blocks in bra/ket sides
            
            auto nbra = braGtoPairsContainer->getNumberOfGtoPairsBlocks();
            
            auto nket = ketGtoPairsContainer->getNumberOfGtoPairsBlocks();
            
            // determine symmetry of bra/ket GTOs pairs containers
            
            auto symbk = ((*braGtoPairsContainer) == (*ketGtoPairsContainer));
            
            // screeners counter
            
            int32_t idx = 0;
            
            // loop over pairs of GTOs blocks
            
            for (int32_t i = (nbra - 1); i >= 0; i--)
            {
                auto bpairs = braGtoPairsContainer->getGtoPairsBlock(i);
                
                auto joff = (symbk) ? i : 0;
                
                for (int32_t j = (nket - 1); j >= joff; j--)
                {
                    #pragma omp task firstprivate(j, idx)
                    {
                        auto kpairs = ketGtoPairsContainer->getGtoPairsBlock(j);
                        
                        auto qqdat = screeningContainer->getScreener(idx);
                        
                        CTwoIntsDistribution distpat(pfock, pden);
                        
                        distpat.setFockContainer(bpairs, kpairs);
                        
                        _compElectronRepulsionForGtoPairsBlocks(distpat, qqdat,
                                                                bpairs, kpairs);
                        
                        // accumulate AO Fock matrix
                        
                        #pragma omp critical (fockacc)
                        distpat.accumulate();
                    }
                    
                    // update screeners counter
                    
                    idx++; 
                }
            }
        }
    }
}

void
CElectronRepulsionIntegralsDriver::_compMaxQValuesForGtoPairsBlock(      double*         qValuesBuffer,
                                                                   const CGtoPairsBlock& gtoPairsBlock) const
{
    // determine number of contracted GTOs pairs
    
    auto ngto = gtoPairsBlock.getNumberOfScreenedContrPairs();
    
    // set up empty integrals screener
    
    CCauchySchwarzScreener qqdat;
    
    // loop over contracted GTOs pairs
    
    for (int32_t i = 0; i < ngto; i++)
    {
        auto cpair = gtoPairsBlock.pick(i);
        
        CTwoIntsDistribution cdist(qValuesBuffer, 1, 1, i, dist2e::qvalues);
        
        _compElectronRepulsionForGtoPairsBlocks(cdist, qqdat, cpair, cpair);
    }
}

CVecMemBlock<double>
CElectronRepulsionIntegralsDriver::_getQValuesBuffer(const CGtoPairsContainer& gtoPairsContainer) const
{
    CVecMemBlock<double> buffvec;
    
    auto nppblk = gtoPairsContainer.getNumberOfGtoPairsBlocks();
    
    for (int32_t i = 0; i < nppblk; i++)
    {
        auto cpairs = gtoPairsContainer.getGtoPairsBlock(i);
        
        buffvec.push_back(CMemBlock<double>(cpairs.getNumberOfScreenedContrPairs()));
    }
    
    return buffvec;
}

