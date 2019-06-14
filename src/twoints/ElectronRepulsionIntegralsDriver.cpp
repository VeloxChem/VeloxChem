//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ElectronRepulsionIntegralsDriver.hpp"

#include <sstream>

#include "GtoPairsContainer.hpp"
#include "SphericalMomentum.hpp"
#include "GenFunc.hpp"
#include "TwoIntsFunc.hpp"
#include "AngularMomentum.hpp"
#include "KetHrrFunc.hpp"
#include "BraHrrFunc.hpp"
#include "StringFormat.hpp"
#include "TwoCentersRecursionFunctions.hpp"
#include "GenIntsFunc.hpp"

#include "ElectronRepulsionRecFuncForSX.hpp"
#include "ElectronRepulsionRecFuncForPX.hpp"
#include "ElectronRepulsionRecFuncForDX.hpp"
#include "ElectronRepulsionRecFuncForFF.hpp"
#include "ElectronRepulsionRecFuncForFG.hpp"
#include "ElectronRepulsionRecFuncForGF.hpp"
#include "ElectronRepulsionRecFuncForFH.hpp"
#include "ElectronRepulsionRecFuncForHF.hpp"
#include "ElectronRepulsionRecFuncForFI.hpp"
#include "ElectronRepulsionRecFuncForIF.hpp"
#include "ElectronRepulsionRecFuncForFK.hpp"
#include "ElectronRepulsionRecFuncForKF.hpp"
#include "ElectronRepulsionRecFuncForFL.hpp"
#include "ElectronRepulsionRecFuncForLF.hpp"
#include "ElectronRepulsionRecFuncForGG.hpp"
#include "ElectronRepulsionRecFuncForGH.hpp"
#include "ElectronRepulsionRecFuncForHG.hpp"
#include "ElectronRepulsionRecFuncForGI.hpp"
#include "ElectronRepulsionRecFuncForIG.hpp"
#include "ElectronRepulsionRecFuncForGK.hpp"
#include "ElectronRepulsionRecFuncForKG.hpp"
#include "ElectronRepulsionRecFuncForGL.hpp"
#include "ElectronRepulsionRecFuncForLG.hpp"
#include "ElectronRepulsionRecFuncForHH.hpp"
#include "ElectronRepulsionRecFuncForHI.hpp"
#include "ElectronRepulsionRecFuncForIH.hpp"
#include "ElectronRepulsionRecFuncForHK.hpp"
#include "ElectronRepulsionRecFuncForKH.hpp"
#include "ElectronRepulsionRecFuncForHL.hpp"
#include "ElectronRepulsionRecFuncForLH.hpp"
#include "ElectronRepulsionRecFuncForII.hpp"
#include "ElectronRepulsionRecFuncForIK.hpp"
#include "ElectronRepulsionRecFuncForKI.hpp"
#include "ElectronRepulsionRecFuncForIL.hpp"
#include "ElectronRepulsionRecFuncForLI.hpp"
#include "ElectronRepulsionRecFuncForKK.hpp"
#include "ElectronRepulsionRecFuncForKL.hpp"
#include "ElectronRepulsionRecFuncForLK.hpp"
#include "ElectronRepulsionRecFuncForLL.hpp"

CElectronRepulsionIntegralsDriver::CElectronRepulsionIntegralsDriver(MPI_Comm comm)
{
    _locRank  = mpi::rank(comm);
    
    _locNodes = mpi::nodes(comm);

    mpi::duplicate(comm, &_locComm);
}

CElectronRepulsionIntegralsDriver::~CElectronRepulsionIntegralsDriver()
{
    mpi::destroy(&_locComm);
}

void
CElectronRepulsionIntegralsDriver::compute(      CAOFockMatrix&       aoFockMatrix,
                                           const CAODensityMatrix&    aoDensityMatrix,
                                           const CMolecule&           molecule,
                                           const CMolecularBasis&     aoBasis,
                                           const CScreeningContainer& screeningContainer) const
{
    // generate GTOs pairs blocks for AO basis on bra side
    
    CGtoPairsContainer bgtopairs(molecule, aoBasis, 1.0e-15);
    
    // split GTOs pairs into batches on bra side
    
    auto bbpairs = bgtopairs.split(_locNodes);
    
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
                                           const CMolecularBasis& aoBasis) const
{
    // generate GTOs pairs blocks for AO basis on bra side
    
    CGtoPairsContainer bgtopairs(molecule, aoBasis, 1.0e-15);
    
    // split GTOs pairs into batches on bra side
    
    auto bbpairs = bgtopairs.split(_locNodes);
    
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
    
    // generate vertical recursion map
    
    auto vrrmap = _setVerticalRecursionMap(tpvec, pmax);

    auto nblk = vrrmap.getNumberOfComponents();
    
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
    
    auto bord = brapairs.getBraAngularMomentum() + brapairs.getKetAngularMomentum()
              + ketpairs.getBraAngularMomentum() + ketpairs.getKetAngularMomentum();
    
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
            
            _compPrimElectronRepulsionInts(pbuffer, vrrmap, bftab,
                                           bargs, bvals, bord, rfacts, rpq, rwp,
                                           rwq, brapairs, ddpairs, nqpdim, i);
            
            // contract primitive electron repulsion integrals
            
            genfunc::contract(khrrbuffer, pbuffer, khrrvec, khrridx, vrrmap, brapairs, ddpairs, nqpdim, nqcdim, i);
            
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
        
            _compPrimElectronRepulsionInts(pbuffer, vrrmap, bftab,
                                           bargs, bvals, bord, rfacts, rpq, rwp,
                                           rwq, brapairs, qqpairs, nqpdim, i);
        
            // contract primitive electron repulsion integrals
        
            genfunc::contract(khrrbuffer, pbuffer, khrrvec, khrridx, vrrmap, brapairs, qqpairs, nqpdim, nqcdim, i);
        
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

CRecursionMap
CElectronRepulsionIntegralsDriver::_setVerticalRecursionMap(const CVecThreeIndexes& leadTerms,
                                                            const int32_t           maxNumberOfPrimPairs) const
{
    CRecursionFunctionsList recfuncs;
    
    recfuncs.add(CRecursionFunction({"Electron Repulsion"}, &t2crecfunc::obRecursionForElectronRepulsion));
    
    CRecursionMap recmap(recblock::cc, maxNumberOfPrimPairs);
    
    for (size_t i = 0; i < leadTerms.size(); i++)
    {
        auto rterm = gintsfunc::genIntegral({"Electron Repulsion"}, leadTerms[i].first(), leadTerms[i].second(), leadTerms[i].third());
        
        recmap.append(gintsfunc::genRecursionMap(rterm, recblock::cc, maxNumberOfPrimPairs, recfuncs));
    }
    
    return recmap;
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
                                                                  const CRecursionMap&        recursionMap, 
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
    erirecfunc::compElectronRepulsionForSSSS(primBuffer, recursionMap, bfTable, bfArguments, bfValues, bfOrder, osFactors, pqDistances,
                                             braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSSSP(primBuffer, recursionMap, wqDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSPSS(primBuffer, recursionMap, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSSSD(primBuffer, recursionMap, osFactors, wqDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSDSS(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSSSF(primBuffer, recursionMap, osFactors, wqDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSFSS(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSSSG(primBuffer, recursionMap, osFactors, wqDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSGSS(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSSSH(primBuffer, recursionMap, osFactors, wqDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSHSS(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSSSI(primBuffer, recursionMap, osFactors, wqDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSISS(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSSSK(primBuffer, recursionMap, osFactors, wqDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSKSS(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSSSL(primBuffer, recursionMap, osFactors, wqDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSLSS(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSPSP(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSPSD(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSDSP(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSPSF(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSFSP(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSPSG(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSGSP(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSPSH(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSHSP(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSPSI(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSISP(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSPSK(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSKSP(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSPSL(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSLSP(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSDSD(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSDSF(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSFSD(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSDSG(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSGSD(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSDSH(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSHSD(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSDSI(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSISD(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSDSK(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSKSD(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSDSL(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSLSD(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSFSF(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSFSG(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSGSF(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSFSH(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSHSF(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSFSI(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSISF(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSFSK(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSKSF(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSFSL(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSLSF(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSGSG(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSGSH(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSHSG(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSGSI(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSISG(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSGSK(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSKSG(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSGSL(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSLSG(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSHSH(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSHSI(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSISH(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSHSK(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSKSH(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSHSL(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSLSH(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSISI(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSISK(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSKSI(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSISL(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSLSI(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSKSK(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSKSL(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSLSK(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    erirecfunc::compElectronRepulsionForSLSL(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
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
CElectronRepulsionIntegralsDriver::_compElectronRepulsionIntegrals(      CAOFockMatrix&       aoFockMatrix,
                                                                   const CAODensityMatrix&    aoDensityMatrix,
                                                                   const CGtoPairsContainer*  braGtoPairsContainer,
                                                                   const CGtoPairsContainer*  ketGtoPairsContainer,
                                                                   const CScreeningContainer* screeningContainer) const
{
    // set up pointers to AO Fock and AO density matrices
    
    auto pfock = &aoFockMatrix;
    
    auto pden = &aoDensityMatrix;
    
    // determine number of GTOs pairs blocks in bra/ket sides
    
    auto nbra = braGtoPairsContainer->getNumberOfGtoPairsBlocks();
    
    auto nket = ketGtoPairsContainer->getNumberOfGtoPairsBlocks();
    
    // determine symmetry of bra/ket GTOs pairs containers
    
    auto symbk = ((*braGtoPairsContainer) == (*ketGtoPairsContainer));
    
    // initialize tasks grid for each MPI process
    
    auto nodpatt = _setTasksGrid(nbra, nket, symbk);
    
    auto ptgrid = nodpatt.data();
    
    #pragma omp parallel shared(braGtoPairsContainer, ketGtoPairsContainer, pfock, pden, nbra, nket, symbk, ptgrid)
    {
        #pragma omp single nowait
        {
            // screeners counter
            
            int32_t idx = 0;
            
            // loop over pairs of GTOs blocks
            
            for (int32_t i = (nbra - 1); i >= 0; i--)
            {
                auto bpairs = braGtoPairsContainer->getGtoPairsBlock(i);
                
                auto joff = (symbk) ? i : 0;
                
                for (int32_t j = (nket - 1); j >= joff; j--)
                {
                    if (ptgrid[idx] == 1)
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

CMemBlock<int32_t>
CElectronRepulsionIntegralsDriver::_setTasksGrid(const int32_t nBraGtoPairsBlocks,
                                                 const int32_t nKetGtoPairsBlocks,
                                                 const bool    isBraEqualKet) const
{
    // determine number of tasks
    
    int32_t nelem = 0;
    
    if (isBraEqualKet)
    {
        nelem = nBraGtoPairsBlocks * (nBraGtoPairsBlocks + 1) / 2;
    }
    else
    {
        nelem = nBraGtoPairsBlocks * nKetGtoPairsBlocks;
    }
    
    // initialize tasks grid
    
    CMemBlock<int32_t> tgrid(nelem);
    
    if (_locNodes == 1)
    {
        mathfunc::set_to(tgrid.data(), 1, nelem);
    }
    else
    {
        mathfunc::set_to(tgrid.data(), 0, nelem);
        
        auto ndim = nelem / _locNodes;
        
        for (int32_t i = 0; i < ndim; i++)
        {
            tgrid.at(i * _locNodes + _locRank) = 1;
        }
        
        auto nrem = nelem % _locNodes;
        
        if (_locRank < nrem)
        {
            tgrid.at(ndim * _locNodes + _locRank) = 1;
        }
    }
    
    return tgrid;
}
