//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ElectricDipoleIntegralsDriver.hpp"

#include "SystemClock.hpp"
#include "AngularMomentum.hpp"
#include "GenFunc.hpp"
#include "MemBlock.hpp"
#include "OneIntsFunc.hpp"
#include "OneIntsDistributor.hpp"

CElectricDipoleIntegralsDriver::CElectricDipoleIntegralsDriver(const int32_t  globRank,
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

CElectricDipoleIntegralsDriver::~CElectricDipoleIntegralsDriver()
{
    
}

CElectricDipoleMatrix
CElectricDipoleIntegralsDriver::compute(const CMolecule&       molecule,
                                        const CMolecularBasis& basis,
                                              MPI_Comm         comm) const
{
    CSystemClock timer;
    
    CElectricDipoleMatrix dipmat;
    
    if (_locRank == mpi::master())
    {
        // set up GTOs container
        
        CGtoContainer bracontr(molecule, basis);
        
        // compute electric dipole integrals
        
        dipmat = _compElectricDipoleIntegrals(&bracontr, &bracontr);
    }
    
    return dipmat;
}

CElectricDipoleMatrix
CElectricDipoleIntegralsDriver::compute(const CMolecule&       molecule,
                                        const CMolecularBasis& braBasis,
                                        const CMolecularBasis& ketBasis,
                                              MPI_Comm         comm) const
{
    CElectricDipoleMatrix dipmat;
    
    if (_locRank == mpi::master())
    {
        // set up GTOs containers
        
        CGtoContainer bracontr(molecule, braBasis);
        
        CGtoContainer ketcontr(molecule, ketBasis);
        
        // compute electric dipole integrals
        
        dipmat = _compElectricDipoleIntegrals(&bracontr, &ketcontr);
    }
    
    return dipmat;
}

CElectricDipoleMatrix
CElectricDipoleIntegralsDriver::compute(const CMolecule&       braMolecule,
                                        const CMolecule&       ketMolecule,
                                        const CMolecularBasis& basis,
                                              MPI_Comm         comm) const
{
    CElectricDipoleMatrix dipmat;
    
    if (_locRank == mpi::master())
    {
        // set up GTOs containers
        
        CGtoContainer bracontr(braMolecule, basis);
        
        CGtoContainer ketcontr(ketMolecule, basis);
        
        // compute electric dipole integrals
        
        dipmat = _compElectricDipoleIntegrals(&bracontr, &ketcontr);
    }
    
    return dipmat;
}

CElectricDipoleMatrix
CElectricDipoleIntegralsDriver::compute(const CMolecule&       braMolecule,
                                        const CMolecule&       ketMolecule,
                                        const CMolecularBasis& braBasis,
                                        const CMolecularBasis& ketBasis,
                                              MPI_Comm         comm) const
{
    CElectricDipoleMatrix dipmat;
    
    if (_locRank == mpi::master())
    {
        // set up GTOs containers
        
        CGtoContainer bracontr(braMolecule, braBasis);
        
        CGtoContainer ketcontr(ketMolecule, ketBasis);
        
        // compute electric integrals
        
        dipmat = _compElectricDipoleIntegrals(&bracontr, &ketcontr);
    }
    
    return dipmat;
}

void
CElectricDipoleIntegralsDriver::compute(      double*    intsBatchX,
                                              double*    intsBatchY,
                                              double*    intsBatchZ,
                                        const CGtoBlock& braGtoBlock,
                                        const CGtoBlock& ketGtoBlock) const
{
    // determine dimensions of integrals batch
    
    auto nrow = angmom::to_SphericalComponents(braGtoBlock.getAngularMomentum())
    
              * braGtoBlock.getNumberOfContrGtos();
    
    auto ncol = angmom::to_SphericalComponents(ketGtoBlock.getAngularMomentum())
    
              * ketGtoBlock.getNumberOfContrGtos();
    
    // set up distribution pattern
    
    COneIntsDistribution distx(intsBatchX, nrow, ncol, dist1e::batch);
    
    COneIntsDistribution disty(intsBatchY, nrow, ncol, dist1e::batch);
    
    COneIntsDistribution distz(intsBatchZ, nrow, ncol, dist1e::batch);
    
    // compute electric dipole integrals
    
    _compElectricDipoleForGtoBlocks(&distx, &disty, &distz, braGtoBlock,
                                    ketGtoBlock);
}

CElectricDipoleMatrix
CElectricDipoleIntegralsDriver::_compElectricDipoleIntegrals(const CGtoContainer* braGtoContainer,
                                                             const CGtoContainer* ketGtoContainer) const
{
    // check if GTOs containers are same on bra and ket sides
    
    auto symbk = ((*braGtoContainer) == (*ketGtoContainer));
    
    // determine dimensions of overlap matrix
    
    auto nrow = braGtoContainer->getNumberOfAtomicOrbitals();
    
    auto ncol = ketGtoContainer->getNumberOfAtomicOrbitals();
    
    // allocate dense matrix for overlap integrals
    
    CDenseMatrix dipxmat(nrow, ncol);
    
    CDenseMatrix dipymat(nrow, ncol);
    
    CDenseMatrix dipzmat(nrow, ncol);
    
    // set up distributio pattern
    
    dist1e dstyp = (symbk) ? dist1e::symsq : dist1e::rect;
    
    COneIntsDistribution* distpatx = new COneIntsDistribution(dipxmat.values(),
                                                              nrow, ncol, dstyp);
    
    COneIntsDistribution* distpaty = new COneIntsDistribution(dipxmat.values(),
                                                              nrow, ncol, dstyp);
    
    COneIntsDistribution* distpatz = new COneIntsDistribution(dipxmat.values(),
                                                              nrow, ncol, dstyp);
    
    // compute overlap integral blocks
    
    #pragma omp parallel shared(braGtoContainer, ketGtoContainer, distpatx, distpaty, distpatz, symbk)
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
                        
                        _compElectricDipoleForGtoBlocks(distpatx, distpaty,
                                                        distpatz, bgtos, kgtos);
                    }
                }
            }
        }
    }
    
    // deallocate distribution pattern
    
    delete distpatx;
    
    delete distpaty;
    
    delete distpatz;
    
    return CElectricDipoleMatrix(dipxmat, dipymat, dipzmat, 0.0, 0.0, 0.0);
}

void
CElectricDipoleIntegralsDriver::_compElectricDipoleForGtoBlocks(      COneIntsDistribution* distPatternX,
                                                                      COneIntsDistribution* distPatternY,
                                                                      COneIntsDistribution* distPatternZ,
                                                                const CGtoBlock&            braGtoBlock,
                                                                const CGtoBlock&            ketGtoBlock) const
{
    // copy GTOs blocks for bra and ket sides
    
    auto bragtos = braGtoBlock;
    
    auto ketgtos = ketGtoBlock;
    
    // copy distribution pattern
    
    auto distpatx = *distPatternX;
    
    auto distpaty = *distPatternY;
    
    auto distpatz = *distPatternZ;
    
    // set up spherical angular momentum for bra and ket sides
    
    CSphericalMomentum bmom(bragtos.getAngularMomentum());
    
    CSphericalMomentum kmom(ketgtos.getAngularMomentum());
    
    // allocate prefactors used in Obara-Saika recursion
    
    auto pdim = ketgtos.getNumberOfPrimGtos();
    
    CMemBlock2D<double> rab(pdim, 3);
    
    auto pmax = bragtos.getMaxContractionDepth();
    
    CMemBlock2D<double> rfacts(pdim, 2 * pmax);
    
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
    
    auto pidx = genfunc::findPairIndex(recidx, recvec, {bang, kang});
    
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
        
        intsfunc::compFactorsForOverlap(rfacts, bragtos, ketgtos, i);
        
        // compute distances: R(PA) = P - A
        
        intsfunc::compDistancesPA(rpa, rab, rfacts, 2, bragtos, ketgtos, i);
        
        // compute distances: R(PB) = P - B
        
        intsfunc::compDistancesPB(rpb, rab, rfacts, 2, bragtos, ketgtos, i);
        
        // compite primitive electric dipole integrals
        
        _compPrimElectricDipoleInts(pbuffer, recvec, recidx, rfacts, rab, rpa, rpb,
                                    bragtos, ketgtos, i);
        
        // contract primitive overlap integrals
        
        genfunc::contract(cartbuffer, pbuffer, pidx, bragtos, ketgtos, i);
        
        // transform Cartesian to spherical integrals
        
        genfunc::transform(spherbuffer, cartbuffer, bmom, kmom, 0, 0, kdim);
        
        // add batch of integrals to integrals matrix
        
        distpatx.distribute(spherbuffer, bragtos, ketgtos, symbk, i);
        
        // add y and z components
    }
}

void
CElectricDipoleIntegralsDriver::_compPrimElectricDipoleInts(      CMemBlock2D<double>&  primBuffer,
                                                            const CVecTwoIndexes&       recPattern,
                                                            const std::vector<int32_t>& recIndexes,
                                                            const CMemBlock2D<double>&  osFactors,
                                                            const CMemBlock2D<double>&  abDistances,
                                                            const CMemBlock2D<double>&  paDistances,
                                                            const CMemBlock2D<double>&  pbDistances,
                                                            const CGtoBlock&            braGtoBlock,
                                                            const CGtoBlock&            ketGtoBlock,
                                                            const int32_t               iContrGto) const
{
   

    // NOTE: add l > 4 recursion here
}

CVecTwoIndexes
CElectricDipoleIntegralsDriver::_getRecursionPattern(const CGtoBlock& braGtoBlock,
                                                     const CGtoBlock& ketGtoBlock) const
{
    // set up angular momentum
    
    auto bang = braGtoBlock.getAngularMomentum();
    
    auto kang = ketGtoBlock.getAngularMomentum();
 
    // set up recursion buffer
    
    CVecTwoIndexes recvec;
    
    recvec.reserve((bang + 1) * (kang + 1));
    
    // set up indexing counters
    
    int32_t spos = 0;
    
    int32_t epos = 1;
    
    // set up initial state of recursion buffer
    
    recvec.push_back(CTwoIndexes(bang, kang));
    
    while (true)
    {
        // internal new recursion terms counter
        
        int32_t nterms = 0;
        
        // generate bra and ket Obara-Saika recursion terms
        
        for (int32_t i = spos; i < epos; i++)
        {
            CTwoIndexes cidx(recvec[i]);
        
            if (cidx.first() != 0)
            {
                // general recursion for bra and ket sides
                
                // (a - 1 | b) term
            
                CTwoIndexes t1idx(cidx.first() - 1,  cidx.second());
            
                if (genfunc::addValidAndUniquePair(recvec, t1idx)) nterms++;
            
                // (a - 2 | b) term
            
                CTwoIndexes t2idx(cidx.first() - 2,  cidx.second());
            
                if (genfunc::addValidAndUniquePair(recvec, t2idx)) nterms++;
            
                // (a - 1 | b - 1) term
            
                CTwoIndexes tkidx(cidx.first() - 1,  cidx.second() - 1);
            
                if (genfunc::addValidAndUniquePair(recvec, tkidx)) nterms++;
            }
            else
            {
                // reduced recursion for ket side
                
                // (0 | b - 1) term
                
                CTwoIndexes t1idx(cidx.first(),  cidx.second() - 1);
                
                if (genfunc::addValidAndUniquePair(recvec, t1idx)) nterms++;
                
                // (0 | b - 2) term
                
                CTwoIndexes t2idx(cidx.first(),  cidx.second() - 2);
                
                if (genfunc::addValidAndUniquePair(recvec, t2idx)) nterms++;
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
CElectricDipoleIntegralsDriver::_getIndexesForRecursionPattern(      std::vector<int32_t>& recIndexes,
                                                               const CVecTwoIndexes&       recPattern,
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
