//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "ThreeCenterElectronRepulsionIntegralsDriver.hpp"

#include <string>

#ifdef MAC_OS_OMP
#include "/opt/intel/compilers_and_libraries/mac/include/omp.h"
#else
#include "omp.h"
#endif

#include "MpiFunc.hpp"
#include "GtoPairsContainer.hpp"
#include "GtoContainer.hpp"
#include "MathFunc.hpp"
#include "StringFormat.hpp"
#include "TwoIntsFunc.hpp"
#include "GenFunc.hpp"
#include "AngularMomentum.hpp"
#include "ThreeCenterEriFunc.hpp"
#include "ThreeCenterHrrFunc.hpp"

CThreeCenterElectronRepulsionIntegralsDriver::CThreeCenterElectronRepulsionIntegralsDriver(const int32_t  globRank,
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

CThreeCenterElectronRepulsionIntegralsDriver::~CThreeCenterElectronRepulsionIntegralsDriver()
{
    
}

void
CThreeCenterElectronRepulsionIntegralsDriver::compute(const CMolecule&       molecule,
                                                      const CMolecularBasis& aoBasis,
                                                      const CMolecularBasis& riBasis,
                                                      const double           threshold, 
                                                            COutputStream&   oStream,
                                                            MPI_Comm         comm) const
{
    CSystemClock eritim;
    
    // generate GTOs pairs blocks for AO basis
    
    CGtoPairsContainer kgtopairs(molecule, aoBasis, 1.0e-13);
    
    // split GTOs pairs into batches
    
    auto kbpairs = kgtopairs.split(5000);
    
    // set up GTOs splitting pattern for RI basis
    
    auto gtopat = _getBatchesOfGtoBlocks(molecule, riBasis, kbpairs);
    
    // generate RI gtos for on each MPI node
    
    CGtoContainer bgtos(molecule, riBasis, gtopat);
    
    // print start header
    
    if (_globRank == mpi::master()) _startHeader(kgtopairs, oStream);
    
    // compute electron repulsion integrals on node
    
    _compElectronRepulsionIntegrals(&bgtos, &kbpairs);
    
    // print evaluation timing statistics
    
    _printTiming(molecule, eritim, oStream);
}

void
CThreeCenterElectronRepulsionIntegralsDriver::compElectronRepulsionForGtoBlocks(const CGtoBlock&      braGtoBlock,
                                                                                const CGtoPairsBlock& ketGtoPairsBlock) const
{
    // copy GTOs and GTOs pairs blocks for bra and ket sides
    
    auto bragtos = braGtoBlock;
    
    auto ketpairs = ketGtoPairsBlock;
    
    // set up spherical angular momentum for bra and ket sides
    
    CSphericalMomentum amom(bragtos.getAngularMomentum());
    
    CSphericalMomentum cmom(ketpairs.getBraAngularMomentum());
    
    CSphericalMomentum dmom(ketpairs.getKetAngularMomentum());
    
    // allocate prefactors used in Obara-Saika recursion
    
    auto pdim = ketpairs.getNumberOfScreenedPrimPairs();
    
    CMemBlock2D<double> raq(pdim, 3);
    
    auto pmax = bragtos.getMaxContractionDepth();
    
    CMemBlock2D<double> rfacts(pdim, 5 * pmax);
    
    CMemBlock2D<double> rw(pdim, 3 * pmax);
    
    CMemBlock2D<double> rwa(pdim, 3 * pmax);
    
    CMemBlock2D<double> rwq(pdim, 3 * pmax);
    
    // generate horizontal recursion pattern
    
    auto hrrvec = _getHorizontalRecursionPattern(bragtos, ketpairs);
    
    // generate intermediate integrals list
    
    auto t0vec = genfunc::getPairsFromTripleIndexes(hrrvec);
    
    // generate vertical recursion pattern
    
    auto vrrvec = _getVerticalRecursionPattern(t0vec);
    
    // set up primitives buffer indexes
    
    std::vector<int32_t> vrridx;
    
    auto nblk = _getIndexesForVerticalRecursionPattern(vrridx, vrrvec, pmax);
    
    // allocate primitives integrals buffer
    
    CMemBlock2D<double> pbuffer(pdim, nblk);
    
    // set up contracted integrals buffer indexes
    
    std::vector<int32_t> t0idx;
    
    nblk = _getIndexesForContractedIntegrals(t0idx, t0vec);
    
    // allocate contracted integrals buffer
    
    auto cdim = ketpairs.getNumberOfScreenedContrPairs();
    
    CMemBlock2D<double> cbuffer(cdim, nblk);
    
    // set up half-transformed integrals buffer indexes
    
    std::vector<int32_t> hrridx;
    
    nblk = _getIndexesForHalfTransformedIntegrals(hrridx, hrrvec);
    
    CMemBlock2D<double> hrrbuffer(cdim, nblk);
    
    // set up angular momentum for bra and ket sides
    
    auto baang = bragtos.getAngularMomentum();
    
    auto kcang = ketpairs.getBraAngularMomentum();
    
    auto kdang = ketpairs.getKetAngularMomentum();
    
    // allocate spherical integrals buffer
    
    auto bcomp = angmom::to_SphericalComponents(bragtos.getAngularMomentum());
    
    nblk = bcomp * angmom::to_SphericalComponents(kcang, kdang);
    
    CMemBlock2D<double> spherbuffer(cdim, nblk);
    
    // set up half transformed integrals position
    
    auto cidx = genfunc::findTripleIndex(hrridx, hrrvec, {baang, kcang, kdang});
    
    // initialize R(CD) = C - D distance for horizontal recursion
    
    auto rcd = ketpairs.getDistancesAB();
    
    // initialize Boys function evaluator
    
    auto bord = genfunc::maxOrderOfPair(vrrvec, 0, 0);
    
    CBoysFunction bftab(bord);
    
    CMemBlock<double> bargs(pdim);
    
    CMemBlock2D<double> bvals(pdim, bord + 1);
    
    // loop over contracted GTOs ob bra side
    
    for (int32_t i = 0; i < bragtos.getNumberOfContrGtos(); i++)
    {
        // compute distances: R(AQ) = A - Q
        
        twointsfunc::compDistancesAQ(raq, bragtos, ketpairs, i);
        
        // compute Obara-Saika recursion factors
        
        twointsfunc::compFactorsForThreeCenterElectronRepulsion(rfacts, bragtos,
                                                                ketpairs, i);
        
        // compute coordinates of center W
        
        twointsfunc::compCoordinatesForW(rw, rfacts, 5, bragtos, ketpairs, i);
        
        // compute distances: R(WA) = W - A
        
        twointsfunc::compDistancesWA(rwa, rw, bragtos, ketpairs, i);
        
        // compute distances: R(WQ) = W - Q;
        
        twointsfunc::compDistancesWQ(rwq, rw, bragtos, ketpairs, i);
        
        // compute primitive electron repulsion integrals
        
        _compPrimElectronRepulsionInts(pbuffer, vrrvec, vrridx, bftab, bargs,
                                       bvals, bord, rfacts, raq, rwa, rwq,
                                       bragtos, ketpairs, i);
        
        // contract primitive electron repulsion integrals
        
        genfunc::contract(cbuffer, pbuffer, t0vec, t0idx, vrrvec, vrridx,
                          bragtos, ketpairs, i);
        
        // transform bra side to spherical form
        
        genfunc::transform(hrrbuffer, cbuffer, amom, hrrvec, hrridx, t0vec,
                           t0idx, cdim);
        
        // apply horizontal recursion
        
        _compContrElectronRepulsionInts(hrrbuffer, hrrvec, hrridx, rcd, bragtos,
                                        ketpairs);
        
        // transform ket side to spherical form
        
        genfunc::transform(spherbuffer, hrrbuffer, cmom, dmom, 0, cidx, cdim, bcomp);
        
        // FIX ME:  distribute or store integrals
    }
}

void
CThreeCenterElectronRepulsionIntegralsDriver::_compElectronRepulsionIntegrals(const CGtoContainer*      braGtoContainer,
                                                                              const CGtoPairsContainer* ketGtoPairsContainer) const
{
    #pragma omp parallel shared(braGtoContainer, ketGtoPairsContainer)
    {
        #pragma omp single nowait
        {
            // determine number of GTOs and GTOs pairs blocks in bra/ket sides
            
            auto nbra = braGtoContainer->getNumberOfGtoBlocks();
            
            auto nket = ketGtoPairsContainer->getNumberOfGtoPairsBlocks();
            
            // loop over pairs of GTOs blocks
            
            for (int32_t i = 0; i < nbra; i++)
            {
                auto bgtos = braGtoContainer->getGtoBlock(i);
                
                for (int32_t j = 0; j < nket; j++)
                {
                    #pragma omp task firstprivate(j)
                    {
                        auto kpairs = ketGtoPairsContainer->getGtoPairsBlock(j);
                        
                        compElectronRepulsionForGtoBlocks(bgtos, kpairs);
                    }
                }
            }
        }
    }
}


void
CThreeCenterElectronRepulsionIntegralsDriver::_compPrimElectronRepulsionInts(      CMemBlock2D<double>&  primBuffer,
                                                                             const CVecThreeIndexes&     recPattern,
                                                                             const std::vector<int32_t>& recIndexes,
                                                                             const CBoysFunction&        bfTable,
                                                                                   CMemBlock<double>&    bfArguments,
                                                                                   CMemBlock2D<double>&  bfValues,
                                                                             const int32_t               bfOrder,
                                                                             const CMemBlock2D<double>&  osFactors,
                                                                             const CMemBlock2D<double>&  aqDistances,
                                                                             const CMemBlock2D<double>&  waDistances,
                                                                             const CMemBlock2D<double>&  wqDistances,
                                                                             const CGtoBlock&            braGtoBlock,
                                                                             const CGtoPairsBlock&       ketGtoPairsBlock,
                                                                             const int32_t               iContrGto) const
{
    // compute (s|g(r,r')|ss) integrals
    
    t3erifunc::compElectronRepulsionForSSS(primBuffer, recPattern, recIndexes,
                                           bfTable, bfArguments, bfValues, bfOrder,
                                           osFactors, aqDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
    
    // compute (s|g(r,r')|sp) integrals
    
    t3erifunc::compElectronRepulsionForSSP(primBuffer, recPattern, recIndexes,
                                           wqDistances, braGtoBlock, ketGtoPairsBlock,
                                           iContrGto);
    
    // compute (p|g(r,r')|ss) integrals
    
    t3erifunc::compElectronRepulsionForPSS(primBuffer, recPattern, recIndexes,
                                           waDistances, braGtoBlock, ketGtoPairsBlock,
                                           iContrGto);
    
    // compute (p|g(r,r')|sp) integrals
    
    t3erifunc::compElectronRepulsionForPSP(primBuffer, recPattern, recIndexes,
                                           osFactors, waDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
    
    // compute (s|g(r,r')|sd) integrals
    
    t3erifunc::compElectronRepulsionForSSD(primBuffer, recPattern, recIndexes,
                                           osFactors, wqDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
    
    // compute (d|g(r,r')|ss) integrals
    
    t3erifunc::compElectronRepulsionForDSS(primBuffer, recPattern, recIndexes,
                                           osFactors, waDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
    
    // compute (p|g(r,r')|sd) integrals
    
    t3erifunc::compElectronRepulsionForPSD(primBuffer, recPattern, recIndexes,
                                           osFactors, waDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
    
    // compute (d|g(r,r')|sp) integrals
    
    t3erifunc::compElectronRepulsionForDSP(primBuffer, recPattern, recIndexes,
                                           osFactors, waDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
    
    // compute (d|g(r,r')|sd) integrals
    
    t3erifunc::compElectronRepulsionForDSD(primBuffer, recPattern, recIndexes,
                                           osFactors, waDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
    
    // compute (s|g(r,r')|sf) integrals
    
    t3erifunc::compElectronRepulsionForSSF(primBuffer, recPattern, recIndexes,
                                           osFactors, wqDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
    
    // compute (f|g(r,r')|ss) integrals
    
    t3erifunc::compElectronRepulsionForFSS(primBuffer, recPattern, recIndexes,
                                           osFactors, waDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
    
    // compute (p|g(r,r')|sf) integrals
    
    t3erifunc::compElectronRepulsionForPSF(primBuffer, recPattern, recIndexes,
                                           osFactors, waDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
    
    // compute (f|g(r,r')|sp) integrals
    
    t3erifunc::compElectronRepulsionForFSP(primBuffer, recPattern, recIndexes,
                                           osFactors, waDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
    
    // compute (d|g(r,r')|sf) integrals
    
    t3erifunc::compElectronRepulsionForDSF(primBuffer, recPattern, recIndexes,
                                           osFactors, waDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
    
    // compute (f|g(r,r')|sd) integrals
    
    t3erifunc::compElectronRepulsionForFSD(primBuffer, recPattern, recIndexes,
                                           osFactors, waDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
    
    // compute (f|g(r,r')|sf) integrals
    
    t3erifunc::compElectronRepulsionForFSF(primBuffer, recPattern, recIndexes,
                                           osFactors, waDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
    
    // compute (s|g(r,r')|sg) integrals
    
    t3erifunc::compElectronRepulsionForSSG(primBuffer, recPattern, recIndexes,
                                           osFactors, wqDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
    
    // compute (g|g(r,r')|ss) integrals
    
    t3erifunc::compElectronRepulsionForGSS(primBuffer, recPattern, recIndexes,
                                           osFactors, waDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
    
    // compute (p|g(r,r')|sg) integrals
    
    t3erifunc::compElectronRepulsionForPSG(primBuffer, recPattern, recIndexes,
                                           osFactors, waDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
    
    // compute (g|g(r,r')|sp) integrals
    
    t3erifunc::compElectronRepulsionForGSP(primBuffer, recPattern, recIndexes,
                                           osFactors, waDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
    
    // compute (d|g(r,r')|sg) integrals
    
    t3erifunc::compElectronRepulsionForDSG(primBuffer, recPattern, recIndexes,
                                           osFactors, waDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
    
    // compute (g|g(r,r')|sd) integrals
    
    t3erifunc::compElectronRepulsionForGSD(primBuffer, recPattern, recIndexes,
                                           osFactors, waDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
    
    // compute (f|g(r,r')|sg) integrals
    
    t3erifunc::compElectronRepulsionForFSG(primBuffer, recPattern, recIndexes,
                                           osFactors, waDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
    
    // compute (g|g(r,r')|sf) integrals
    
    t3erifunc::compElectronRepulsionForGSF(primBuffer, recPattern, recIndexes,
                                           osFactors, waDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
    
    // compute (g|g(r,r')|sg) integrals
    
    t3erifunc::compElectronRepulsionForGSG(primBuffer, recPattern, recIndexes,
                                           osFactors, waDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
    
    // compute (s|g(r,r')|sh) integrals
    
    t3erifunc::compElectronRepulsionForSSH(primBuffer, recPattern, recIndexes,
                                           osFactors, wqDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
    
    // compute (p|g(r,r')|sh) integrals
    
    t3erifunc::compElectronRepulsionForPSH(primBuffer, recPattern, recIndexes,
                                           osFactors, waDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
    
    // compute (d|g(r,r')|sh) integrals
    
    t3erifunc::compElectronRepulsionForDSH(primBuffer, recPattern, recIndexes,
                                           osFactors, waDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
    
    // compute (f|g(r,r')|sh) integrals
    
    t3erifunc::compElectronRepulsionForFSH(primBuffer, recPattern, recIndexes,
                                           osFactors, waDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
    
    // compute (g|g(r,r')|sh) integrals
    
    t3erifunc::compElectronRepulsionForGSH(primBuffer, recPattern, recIndexes,
                                           osFactors, waDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
    
    // compute (s|g(r,r')|si) integrals
    
    t3erifunc::compElectronRepulsionForSSI(primBuffer, recPattern, recIndexes,
                                           osFactors, wqDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
    
    // compute (p|g(r,r')|si) integrals
    
    t3erifunc::compElectronRepulsionForPSI(primBuffer, recPattern, recIndexes,
                                           osFactors, waDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
    
    // compute (d|g(r,r')|si) integrals
    
    t3erifunc::compElectronRepulsionForDSI(primBuffer, recPattern, recIndexes,
                                           osFactors, waDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
    
    // compute (f|g(r,r')|si) integrals
    
    t3erifunc::compElectronRepulsionForFSI(primBuffer, recPattern, recIndexes,
                                           osFactors, waDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
    
    // compute (g|g(r,r')|si) integrals
    
    t3erifunc::compElectronRepulsionForGSI(primBuffer, recPattern, recIndexes,
                                           osFactors, waDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
    
    // compute (s|g(r,r')|sk) integrals
    
    t3erifunc::compElectronRepulsionForSSK(primBuffer, recPattern, recIndexes,
                                           osFactors, wqDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
    
    // compute (p|g(r,r')|sk) integrals
    
    t3erifunc::compElectronRepulsionForPSK(primBuffer, recPattern, recIndexes,
                                           osFactors, waDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
    
    // compute (d|g(r,r')|sk) integrals
    
    t3erifunc::compElectronRepulsionForDSK(primBuffer, recPattern, recIndexes,
                                           osFactors, waDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
    
    // compute (f|g(r,r')|sk) integrals
    
    t3erifunc::compElectronRepulsionForFSK(primBuffer, recPattern, recIndexes,
                                           osFactors, waDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
    
    // compute (g|g(r,r')|sk) integrals
    
    t3erifunc::compElectronRepulsionForGSK(primBuffer, recPattern, recIndexes,
                                           osFactors, waDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
    
    // compute (s|g(r,r')|sl) integrals
    
    t3erifunc::compElectronRepulsionForSSL(primBuffer, recPattern, recIndexes,
                                           osFactors, wqDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
    
    // compute (p|g(r,r')|sl) integrals
    
    t3erifunc::compElectronRepulsionForPSL(primBuffer, recPattern, recIndexes,
                                           osFactors, waDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
    
    // compute (d|g(r,r')|sl) integrals
    
    t3erifunc::compElectronRepulsionForDSL(primBuffer, recPattern, recIndexes,
                                           osFactors, waDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
    
    // compute (f|g(r,r')|sl) integrals
    
    t3erifunc::compElectronRepulsionForFSL(primBuffer, recPattern, recIndexes,
                                           osFactors, waDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
    
    // compute (g|g(r,r')|sl) integrals
    
    t3erifunc::compElectronRepulsionForGSL(primBuffer, recPattern, recIndexes,
                                           osFactors, waDistances, braGtoBlock,
                                           ketGtoPairsBlock, iContrGto);
}

void
CThreeCenterElectronRepulsionIntegralsDriver::_compContrElectronRepulsionInts(      CMemBlock2D<double>&  contrBuffer,
                                                                              const CVecThreeIndexes&     recPattern,
                                                                              const std::vector<int32_t>& recIndexes,
                                                                              const CMemBlock2D<double>&  cdDistances,
                                                                              const CGtoBlock&            braGtoBlock,
                                                                              const CGtoPairsBlock&       ketGtoPairsBlock) const
{
    // set up angular momentum on bra side
    
    auto bang = braGtoBlock.getAngularMomentum();
    
    // compute (x|g(r,r')|pp) integrals

    t3hrrfunc::compElectronRepulsionForXPP(contrBuffer, recPattern, recIndexes,
                                           cdDistances, bang, ketGtoPairsBlock);
    
    // compute (x|g(r,r')|pd) integrals
    
    t3hrrfunc::compElectronRepulsionForXPD(contrBuffer, recPattern, recIndexes,
                                           cdDistances, bang, ketGtoPairsBlock);
    
    // compute (x|g(r,r')|pf) integrals
    
    t3hrrfunc::compElectronRepulsionForXPF(contrBuffer, recPattern, recIndexes,
                                           cdDistances, bang, ketGtoPairsBlock);
    
    // compute (x|g(r,r')|pg) integrals
    
    t3hrrfunc::compElectronRepulsionForXPG(contrBuffer, recPattern, recIndexes,
                                           cdDistances, bang, ketGtoPairsBlock);
    
    // compute (x|g(r,r')|ph) integrals
    
    t3hrrfunc::compElectronRepulsionForXPH(contrBuffer, recPattern, recIndexes,
                                           cdDistances, bang, ketGtoPairsBlock);
    
    // compute (x|g(r,r')|pi) integrals
    
    t3hrrfunc::compElectronRepulsionForXPI(contrBuffer, recPattern, recIndexes,
                                           cdDistances, bang, ketGtoPairsBlock);
    
    // compute (x|g(r,r')|pk) integrals
    
    t3hrrfunc::compElectronRepulsionForXPK(contrBuffer, recPattern, recIndexes,
                                           cdDistances, bang, ketGtoPairsBlock);
    
    // compute (x|g(r,r')|dd) integrals
    
    t3hrrfunc::compElectronRepulsionForXDD(contrBuffer, recPattern, recIndexes,
                                           cdDistances, bang, ketGtoPairsBlock);
    
    // compute (x|g(r,r')|df) integrals
    
    t3hrrfunc::compElectronRepulsionForXDF(contrBuffer, recPattern, recIndexes,
                                           cdDistances, bang, ketGtoPairsBlock);
    
    // compute (x|g(r,r')|dg) integrals
    
    t3hrrfunc::compElectronRepulsionForXDG(contrBuffer, recPattern, recIndexes,
                                           cdDistances, bang, ketGtoPairsBlock);
    
    // compute (x|g(r,r')|dh) integrals
    
    t3hrrfunc::compElectronRepulsionForXDH(contrBuffer, recPattern, recIndexes,
                                           cdDistances, bang, ketGtoPairsBlock);
    
    // compute (x|g(r,r')|di) integrals
    
    t3hrrfunc::compElectronRepulsionForXDI(contrBuffer, recPattern, recIndexes,
                                           cdDistances, bang, ketGtoPairsBlock);
    
    // compute (x|g(r,r')|ff) integrals
    
    t3hrrfunc::compElectronRepulsionForXFF(contrBuffer, recPattern, recIndexes,
                                           cdDistances, bang, ketGtoPairsBlock);
    
    // compute (x|g(r,r')|fg) integrals
    
    t3hrrfunc::compElectronRepulsionForXFG(contrBuffer, recPattern, recIndexes,
                                           cdDistances, bang, ketGtoPairsBlock);
    
    // compute (x|g(r,r')|fh) integrals
    
    t3hrrfunc::compElectronRepulsionForXFH(contrBuffer, recPattern, recIndexes,
                                           cdDistances, bang, ketGtoPairsBlock);
    
    // compute (x|g(r,r')|gg) integrals
    
    t3hrrfunc::compElectronRepulsionForXGG(contrBuffer, recPattern, recIndexes,
                                           cdDistances, bang, ketGtoPairsBlock);
}

void
CThreeCenterElectronRepulsionIntegralsDriver::_startHeader(const CGtoPairsContainer& gtoPairs,
                                                                 COutputStream&      oStream) const
{
    oStream << fmt::header << "Three-Center Electron Repulsion Integrals" << fmt::end;
    
    oStream << std::string(43, '=') << fmt::end << fmt::blank;
    
    // GTO pairs screening information
    
    gtoPairs.printScreeningInfo(oStream);
    
    oStream << fmt::blank;
}

void
CThreeCenterElectronRepulsionIntegralsDriver::_printTiming(const CMolecule&     molecule,
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
        
        std::string str("Three-Center Integrals Evaluation Timings: ");
        
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

CMemBlock2D<int32_t>
CThreeCenterElectronRepulsionIntegralsDriver::_getBatchesOfGtoBlocks(const CMolecule&          molecule,
                                                                     const CMolecularBasis&    riBasis,
                                                                     const CGtoPairsContainer& gtoPairsContainer) const
{
    // determine dimensions of atoms batch for each MPI node
    
    auto natoms = molecule.getNumberOfAtoms();
    
    auto nodatm = mpi::batch_size(natoms, _locRank, _locNodes);
    
    auto nodoff = mpi::batch_offset(natoms, _locRank, _locNodes);
    
    // determine max. task vector dimensions for single atoms batch
    
    auto mtasks = (riBasis.getMaxAngularMomentum() + 1)
    
                * gtoPairsContainer.getNumberOfGtoPairsBlocks();
    
    // determine max. number of threads
    
    int32_t mthreads = 10 * omp_get_max_threads();
    
    // determine number of atomic batches
    
    auto nblocks = mthreads / mtasks + 1;
    
    CMemBlock2D<int32_t> batches(nblocks, 2);
    
    auto bpos = batches.data(0);
   
    auto bdim = batches.data(1);
    
    for (int32_t i = 0; i < nblocks; i++)
    {
        bdim[i] = mpi::batch_size(nodatm, i, nblocks);
    }
    
    mathfunc::indexes(bpos, bdim, nodoff, nblocks);
    
    return batches;
}

CVecThreeIndexes
CThreeCenterElectronRepulsionIntegralsDriver::_getHorizontalRecursionPattern(const CGtoBlock&      braGtoBlock,
                                                                             const CGtoPairsBlock& ketGtoPairsBlock) const
{
    // set up angular momentum
    
    auto anga = braGtoBlock.getAngularMomentum();
    
    auto angc = ketGtoPairsBlock.getBraAngularMomentum();
    
    auto angd = ketGtoPairsBlock.getKetAngularMomentum();
    
    // set up recursion buffer
    
    CVecThreeIndexes recvec;
    
    recvec.reserve((angc + 1) * (angd + 1));
    
    // set up indexing counters
    
    int32_t spos = 0;
    
    int32_t epos = 1;
    
    // set up initial state of recursion buffer
    
    recvec.push_back(CThreeIndexes(anga, angc, angd));
    
    while (true)
    {
        // internal new recursion terms counter
        
        int32_t nterms = 0;
        
        // generate bra and ket Obara-Saika recursion terms
        
        for (int32_t i = spos; i < epos; i++)
        {
            CThreeIndexes cidx(recvec[i]);
            
            // (a |g(r,r')| (c - 1) d) term
            
            CThreeIndexes t0idx(cidx.first(),  cidx.second() - 1,
                                
                                cidx.third());
            
            if (genfunc::addValidAndUniqueTriple(recvec, t0idx)) nterms++;
            
            // (a |g(r,r')| (c - 1) (d + 1)) term
            
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
CThreeCenterElectronRepulsionIntegralsDriver::_getVerticalRecursionPattern(const CVecThreeIndexes& leadTerms) const
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
                // (a - 1 |g(r,r')| 0 d)^(m+1) term
                
                CThreeIndexes t0idx(cidx.first() - 1,  cidx.second(),
                                    
                                    cidx.third() + 1);
                
                if (genfunc::addValidAndUniqueTriple(recvec, t0idx)) nterms++;
                
                // (a - 2 |g(r,r')| 0 d)^(m) term
                
                CThreeIndexes t1idx(cidx.first() - 2,  cidx.second(),
                                    
                                    cidx.third());
                
                if (genfunc::addValidAndUniqueTriple(recvec, t1idx)) nterms++;
                
                // (a - 2 |g(r,r')| 0 d)^(m+1) term
                
                CThreeIndexes t2idx(cidx.first() - 2,  cidx.second(),
                                    
                                    cidx.third() + 1);
                
                if (genfunc::addValidAndUniqueTriple(recvec, t2idx)) nterms++;
                
                // (a - 1 |g(r,r')| 0 (d - 1))^(m+1) term
                
                CThreeIndexes tkidx(cidx.first() - 1,  cidx.second() - 1,
                                    
                                    cidx.third() + 1);
                
                if (genfunc::addValidAndUniqueTriple(recvec, tkidx)) nterms++;
            }
            else
            {
                // (0 |g(r,r')| 0 (d - 1))^(m) term
                
                CThreeIndexes t0idx(cidx.first(),  cidx.second() - 1,
                                    
                                    cidx.third());
                
                if (genfunc::addValidAndUniqueTriple(recvec, t0idx)) nterms++;
                
                // (0 |g(r,r')| 0 (d - 1))^(m+1) term
                
                CThreeIndexes t1idx(cidx.first(),  cidx.second() - 1,
                                    
                                    cidx.third() + 1);
                
                if (genfunc::addValidAndUniqueTriple(recvec, t1idx)) nterms++;
                
                // (0 |g(r,r')| 0 (d - 2))^(m) term
                
                CThreeIndexes t2idx(cidx.first(),  cidx.second() - 2,
                                    
                                    cidx.third());
                
                if (genfunc::addValidAndUniqueTriple(recvec, t2idx)) nterms++;
                
                // (0 |g(r,r')| 0 (d - 2))^(m+1) term
                
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
CThreeCenterElectronRepulsionIntegralsDriver::_getIndexesForVerticalRecursionPattern(      std::vector<int32_t>& recIndexes,
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

int32_t
CThreeCenterElectronRepulsionIntegralsDriver::_getIndexesForContractedIntegrals(      std::vector<int32_t>& contrIndexes,
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

int32_t
CThreeCenterElectronRepulsionIntegralsDriver::_getIndexesForHalfTransformedIntegrals(std::vector<int32_t>&   intsIndexes,
                                                                                     const CVecThreeIndexes& intsListing) const
{
    // clear vector and reserve memory
    
    intsIndexes.clear();
    
    intsIndexes.reserve(intsListing.size() + 1);
    
    // loop over integrals listing
    
    int32_t nblk = 0;
    
    for (size_t i = 0; i < intsListing.size(); i++)
    {
        intsIndexes.push_back(nblk);
        
        nblk += angmom::to_SphericalComponents(intsListing[i].first())
        
              *  angmom::to_CartesianComponents(intsListing[i].second(),
                                                intsListing[i].third());
    }
    
    return nblk;
}
