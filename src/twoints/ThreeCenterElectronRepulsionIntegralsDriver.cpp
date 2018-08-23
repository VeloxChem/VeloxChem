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
#include "SystemClock.hpp"
#include "GtoPairsContainer.hpp"
#include "GtoContainer.hpp"
#include "MathFunc.hpp"
#include "StringFormat.hpp"
#include "TwoIntsFunc.hpp"
#include "GenFunc.hpp"
#include "AngularMomentum.hpp"
#include "ThreeCenterEriFunc.hpp"

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
    
    // set up GTOs splitting pattern for RI basis
    
    auto gtopat = _getBatchesOfGtoBlocks(molecule, riBasis, kgtopairs);
    
    // generate RI gtos for on each MPI node
    
    CGtoContainer bgtos(molecule, riBasis, gtopat);
    
    // print start header
    
    if (_globRank == mpi::master()) _startHeader(kgtopairs, oStream);
    
    // compute electron repulsion integrals on node
    
    _compElectronRepulsionIntegrals(&bgtos, &kgtopairs);
    
    printf("node: %i ngtos: %i time: %lf sec.\n", _locRank, bgtos.getNumberOfAtomicOrbitals(), eritim.getElapsedTimeInSeconds());
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
                        
                        _compElectronRepulsionForGtoBlocks(bgtos, kpairs);
                    }
                }
            }
        }
    }
}

void
CThreeCenterElectronRepulsionIntegralsDriver::_compElectronRepulsionForGtoBlocks(const CGtoBlock&      braGtoBlock,
                                                                                 const CGtoPairsBlock& ketGtoPairsBlock) const
{
    // copy GTOs and GTOs pairs blocks for bra and ket sides
    
    auto bragtos = braGtoBlock;
    
    auto ketpairs = ketGtoPairsBlock;
 
    // testing print out
    
    printf("*** INTEGRALS (%i|%i,%i):\n", bragtos.getAngularMomentum(),
           ketpairs.getBraAngularMomentum(), ketpairs.getKetAngularMomentum());
    
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
    
    CMemBlock2D<double> rwd(pdim, 3 * pmax);
    
    // generate horizontal recursion pattern
    
    auto hrrvec = _getHorizontalRecursionPattern(bragtos, ketpairs);
    
    // generate intermediate integrals list
    
    auto t0vec = genfunc::getPairsFromTripleIndexes(hrrvec);
    
    // generate vertical recursion pattern
    
    auto vrrvec = _getVerticalRecursionPattern(t0vec, bragtos, ketpairs);
    
    // set up primitives buffer indexes
    
    std::vector<int32_t> vrridx;
    
    auto nblk = _getIndexesForVerticalRecursionPattern(vrridx, vrrvec, pmax);
    
    // allocate primitives integrals buffer
    
    CMemBlock2D<double> pbuffer(pdim, nblk);
    
    // FIX ME: add other buffers
    
    // initialize Boys function evaluator
    
    auto bord = genfunc::maxOrderOfPair(vrrvec, 0, 0);
    
    CBoysFunction bftab(bord);
    
    CMemBlock<double> bargs(pdim);
    
    CMemBlock2D<double> bvals(pdim, bord + 1);
    
    // testing codde
    
    for (int32_t i = 0; i < hrrvec.size(); i++)
    {
        printf("HRR: (%i,%i,%i)\n", hrrvec[i].first(), hrrvec[i].second(),
               hrrvec[i].third());
    }
    
    printf("Intermidiates: ");
    
    for (int32_t i = 0; i < t0vec.size(); i++)
    {
        printf("(%i,%i)^(%i) ", t0vec[i].first(), t0vec[i].second(),
               t0vec[i].third());
    }
    
    printf("\n");
    
    for (int32_t i = 0; i < vrrvec.size(); i++)
    {
        printf("VRR: (%i,%i)^(%i) Index: %i\n", vrrvec[i].first(), vrrvec[i].second(),
               vrrvec[i].third(), vrridx[i]);
    }
    
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
        
        // compute distances: R(WD) = W - D;
        
        twointsfunc::compDistancesWD(rwd, rw, bragtos, ketpairs, i); 
        
        // compute primitive electron repulsion integrals
        
        _compPrimElectronRepulsionInts(pbuffer, vrrvec, vrridx, bftab, bargs,
                                       bvals, bord, rfacts, raq, rwa, rwd,
                                       bragtos, ketpairs, i);
        
        // contract primitive electron repulsion integrals
        
        // transform bra side to spherical form
        
        // apply horizontal recursion
        
        // transform ket side to spherical form
        
        // store computed integrals
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
                                                                             const CMemBlock2D<double>&  wdDistances,
                                                                             const CGtoBlock&            braGtoBlock,
                                                                             const CGtoPairsBlock&       ketGtoPairsBlock,
                                                                             const int32_t               iContrGto) const
{
    // compute (s|g(r,r')|ss) integrals
    
    t3erifunc::compElectronicPotentialForSSS(primBuffer, recPattern, recIndexes,
                                             bfTable, bfArguments, bfValues, bfOrder,
                                             osFactors, aqDistances, braGtoBlock,
                                             ketGtoPairsBlock, iContrGto); 
}

void
CThreeCenterElectronRepulsionIntegralsDriver::_startHeader(const CGtoPairsContainer& gtoPairs,
                                                                 COutputStream&      oStream) const
{
    oStream << fmt::header << "Three-Center Electron Repulsion Integrals" << fmt::end;
    
    oStream << std::string(43, '=') << fmt::end << fmt::blank;
    
    // GTO pairs screening information
    
    gtoPairs.printScreeningInfo(oStream);
    
    // memory storage information
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
CThreeCenterElectronRepulsionIntegralsDriver::_getVerticalRecursionPattern(const CVecThreeIndexes&  leadTerms,
                                                                           const CGtoBlock&         braGtoBlock,
                                                                           const CGtoPairsBlock&    ketGtoPairsBlock) const
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
