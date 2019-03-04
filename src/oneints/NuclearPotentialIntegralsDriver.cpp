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
#include "NuclearPotentialRecFunc.hpp"
#include "StringFormat.hpp"

CNuclearPotentialIntegralsDriver::CNuclearPotentialIntegralsDriver(const int32_t  globRank,
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

CNuclearPotentialIntegralsDriver::~CNuclearPotentialIntegralsDriver()
{
}

CNuclearPotentialMatrix
CNuclearPotentialIntegralsDriver::compute(const CMolecule&       molecule,
                                          const CMolecularBasis& basis,
                                                MPI_Comm         comm) const
{
    CSystemClock timer;
    
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
                                          const CMolecule&       pchgMolecule,
                                                MPI_Comm         comm) const
{
    CSystemClock timer;
    
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
                                          const CMolecule&       pchgMolecule,
                                                MPI_Comm         comm) const
{
    CSystemClock timer;
    
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
                                          const CMolecule&       pchgMolecule,
                                                MPI_Comm         comm) const
{
    CSystemClock timer;
    
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
                                          const CMolecule&       pchgMolecule,
                                                MPI_Comm         comm) const
{
    CSystemClock timer;
    
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
    
    // set up spherical angular momentum for bra and ket sides
    
    CSphericalMomentum bmom(bragtos.getAngularMomentum());
    
    CSphericalMomentum kmom(ketgtos.getAngularMomentum());
    
    // allocate prefactors used in Obara-Saika recursion
    
    auto pdim = ketgtos.getNumberOfPrimGtos();
    
    CMemBlock2D<double> rab(pdim, 3);
    
    auto pmax = bragtos.getMaxContractionDepth();
    
    CMemBlock2D<double> rfacts(pdim, 3 * pmax);
    
    CMemBlock2D<double> rp(pdim, 3 * pmax);
    
    CMemBlock2D<double> rpa(pdim, 3 * pmax);
    
    CMemBlock2D<double> rpb(pdim, 3 * pmax);
    
    CMemBlock2D<double> rpc(pdim, 3 * pmax);
    
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
    
    // allocate accumulation buffer
    
    auto ncart = angmom::to_CartesianComponents(bang, kang);
    
    CMemBlock2D<double> accbuffer(pdim, ncart * pmax);
    
    // set up contracted GTOs dimensions
    
    auto kdim = ketgtos.getNumberOfContrGtos();
    
    // allocate contracted Cartesian integrals buffer
    
    CMemBlock2D<double> cartbuffer(kdim, ncart);
    
    // allocate contracted spherical integrals buffer
    
    auto nspher = angmom::to_SphericalComponents(bang, kang);
    
    CMemBlock2D<double> spherbuffer(kdim, nspher);
    
    // initialize Boys function evaluator
    
    auto bord = genfunc::maxOrderOfPair(recvec, 0, 0);
    
    CBoysFunction bftab(bord);
    
    CMemBlock<double> bargs(pdim);
    
    CMemBlock2D<double> bvals(pdim, bord + 1);
    
    // determine bra and ket sides symmetry
    
    bool symbk = (bragtos == ketgtos);
    
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
        
        // reset accumulation buffer
        
        accbuffer.zero();
        
        // loop over charges
        
        for (int32_t j = 0; j < qvalues.size(); j++)
        {
            // compute distances: R(PC) = P - C
            
            intsfunc::compDistancesPC(rpc, rp, qcoords, bragtos, ketgtos, i, j);
            
            // compute primitive integrals
            
            _compPrimNuclearPotentialInts(pbuffer, recvec, recidx, bftab, bargs,
                                          bvals, bord, rfacts, rab, rpa, rpb,
                                          rpc, bragtos, ketgtos, i);
            
            // add scaled contribution to accumulation buffer
            
            _addPointChargeContribution(accbuffer, pbuffer, pidx, qvalues,
                                        bragtos, ketgtos, i, j);
        }
        
        // contract primitive overlap integrals
        
        genfunc::contract(cartbuffer, accbuffer, 0, bragtos, ketgtos, i);
        
        // transform Cartesian to spherical integrals
        
        genfunc::transform(spherbuffer, cartbuffer, bmom, kmom, 0, 0, kdim);
        
        // add batch of integrals to integrals matrix
        
        // add batch of integrals to integrals matrix
        
        distpat.distribute(spherbuffer, bragtos, ketgtos, symbk, i);
    }
}

void
CNuclearPotentialIntegralsDriver::_compPrimNuclearPotentialInts(      CMemBlock2D<double>&  primBuffer,
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
                                                                const CMemBlock2D<double>&  pcDistances,
                                                                const CGtoBlock&            braGtoBlock,
                                                                const CGtoBlock&            ketGtoBlock,
                                                                const int32_t               iContrGto) const
{
    // compute (s|A(0)|s) integrals
    
    npotrecfunc::compNuclearPotentialForSS(primBuffer, recPattern, recIndexes,
                                           bfTable, bfArguments, bfValues,
                                           bfOrder, osFactors, abDistances,
                                           pcDistances, braGtoBlock, ketGtoBlock,
                                           iContrGto);
    
    // compute (s|A(0)|p) integrals
    
    npotrecfunc::compNuclearPotentialForSP(primBuffer, recPattern, recIndexes,
                                           pbDistances, pcDistances, braGtoBlock,
                                           ketGtoBlock, iContrGto);
    
    // compute (p|A(0)|s) integrals
    
    npotrecfunc::compNuclearPotentialForPS(primBuffer, recPattern, recIndexes,
                                           paDistances, pcDistances, braGtoBlock,
                                           ketGtoBlock, iContrGto);
    
    // compute (p|A(0)|p) integrals
    
    npotrecfunc::compNuclearPotentialForPP(primBuffer, recPattern, recIndexes,
                                           osFactors, paDistances, pcDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto);
    
    // compute (s|A(0)|d) integrals
    
    npotrecfunc::compNuclearPotentialForSD(primBuffer, recPattern, recIndexes,
                                           osFactors, pbDistances, pcDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto);
    
    // compute (d|A(0)|s) integrals
    
    npotrecfunc::compNuclearPotentialForDS(primBuffer, recPattern, recIndexes,
                                           osFactors, paDistances, pcDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto);
    
    // compute (p|A(0)|d) integrals
    
    npotrecfunc::compNuclearPotentialForPD(primBuffer, recPattern, recIndexes,
                                           osFactors, paDistances, pcDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto);
    
    // compute (d|A(0)|p) integrals
    
    npotrecfunc::compNuclearPotentialForDP(primBuffer, recPattern, recIndexes,
                                           osFactors, paDistances, pcDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto);
    
    // compute (d|A(0)|d) integrals
    
    npotrecfunc::compNuclearPotentialForDD(primBuffer, recPattern, recIndexes,
                                           osFactors, paDistances, pcDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto);
    
    // compute (s|A(0)|f) integrals
    
    npotrecfunc::compNuclearPotentialForSF(primBuffer, recPattern, recIndexes,
                                           osFactors, pbDistances, pcDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto);
    
    // compute (f|A(0)|s) integrals
    
    npotrecfunc::compNuclearPotentialForFS(primBuffer, recPattern, recIndexes,
                                           osFactors, paDistances, pcDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto);
    
    // compute (p|A(0)|f) integrals
    
    npotrecfunc::compNuclearPotentialForPF(primBuffer, recPattern, recIndexes,
                                           osFactors, paDistances, pcDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto);
    
    // compute (f|A(0)|p) integrals
    
    npotrecfunc::compNuclearPotentialForFP(primBuffer, recPattern, recIndexes,
                                           osFactors, paDistances, pcDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto);
    
    // compute (d|A(0)|f) integrals
    
    npotrecfunc::compNuclearPotentialForDF(primBuffer, recPattern, recIndexes,
                                           osFactors, paDistances, pcDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto);
    
    // compute (f|A(0)|d) integrals
    
    npotrecfunc::compNuclearPotentialForFD(primBuffer, recPattern, recIndexes,
                                           osFactors, paDistances, pcDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto);
    
    // compute (f|A(0)|f) integrals
    
    npotrecfunc::compNuclearPotentialForFF(primBuffer, recPattern, recIndexes,
                                           osFactors, paDistances, pcDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto);
    
    // compute (s|A(0)|g) integrals
    
    npotrecfunc::compNuclearPotentialForSG(primBuffer, recPattern, recIndexes,
                                           osFactors, pbDistances, pcDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto);
    
    // compute (g|A(0)|s) integrals
    
    npotrecfunc::compNuclearPotentialForGS(primBuffer, recPattern, recIndexes,
                                           osFactors, paDistances, pcDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto);
    
    // compute (p|A(0)|g) integrals
    
    npotrecfunc::compNuclearPotentialForPG(primBuffer, recPattern, recIndexes,
                                           osFactors, paDistances, pcDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto);
    
    // compute (g|A(0)|p) integrals
    
    npotrecfunc::compNuclearPotentialForGP(primBuffer, recPattern, recIndexes,
                                           osFactors, paDistances, pcDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto);
    
    // compute (d|A(0)|g) integrals
    
    npotrecfunc::compNuclearPotentialForDG(primBuffer, recPattern, recIndexes,
                                           osFactors, paDistances, pcDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto);
    
    // compute (g|A(0)|d) integrals
    
    npotrecfunc::compNuclearPotentialForGD(primBuffer, recPattern, recIndexes,
                                           osFactors, paDistances, pcDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto);
    
    // compute (f|A(0)|g) integrals
    
    npotrecfunc::compNuclearPotentialForFG(primBuffer, recPattern, recIndexes,
                                           osFactors, paDistances, pcDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto);
    
    // compute (g|A(0)|f) integrals
    
    npotrecfunc::compNuclearPotentialForGF(primBuffer, recPattern, recIndexes,
                                           osFactors, paDistances, pcDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto);
    
    // compute (g|A(0)|g) integrals
    
    npotrecfunc::compNuclearPotentialForGG(primBuffer, recPattern, recIndexes,
                                           osFactors, paDistances, pcDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto);
    
    // NOTE: add l > 4 recursion here
}

CVecThreeIndexes
CNuclearPotentialIntegralsDriver::_getRecursionPattern(const CGtoBlock& braGtoBlock,
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
            
            // nuclear potential recursion
                
            if (cidx.first() != 0)
            {
                // general recursion for bra and ket sides
                
                // (a - 1 |A(0)| b)^(m) term
                    
                CThreeIndexes t10idx(cidx.first() - 1,  cidx.second(),
                                    
                                     cidx.third());
                    
                if (genfunc::addValidAndUniqueTriple(recvec, t10idx)) nterms++;
                    
                // (a - 1 |A(0)| b)^(m+1) term
                    
                CThreeIndexes t11idx(cidx.first() - 1,  cidx.second(),
                                    
                                     cidx.third() + 1);
                    
                if (genfunc::addValidAndUniqueTriple(recvec, t11idx)) nterms++;
                    
                // (a - 2 |A(0)| b)^(m) term
                
                CThreeIndexes t20idx(cidx.first() - 2,  cidx.second(),
                                     
                                     cidx.third());
                
                if (genfunc::addValidAndUniqueTriple(recvec, t20idx)) nterms++;
                
                // (a - 2 |A(0)| b)^(m+1) term
                
                CThreeIndexes t21idx(cidx.first() - 2,  cidx.second(),
                                     
                                     cidx.third() + 1);
                
                if (genfunc::addValidAndUniqueTriple(recvec, t21idx)) nterms++;
                
                // (a - 1 |A(0)| b - 1)^(m) term
                
                CThreeIndexes tk0idx(cidx.first() - 1,  cidx.second() - 1,
                                    
                                     cidx.third());
                    
                if (genfunc::addValidAndUniqueTriple(recvec, tk0idx)) nterms++;
                
                // (a - 1 |A(0)| b - 1)^(m+1) term
                
                CThreeIndexes tk1idx(cidx.first() - 1,  cidx.second() - 1,
                                     
                                     cidx.third() + 1);
                
                if (genfunc::addValidAndUniqueTriple(recvec, tk1idx)) nterms++;
                }
                else
                {
                    // simplified recursion for ket sides
                    
                    // (0 |A(0)| b - 1)^(m) term
                    
                    CThreeIndexes t10idx(cidx.first(),  cidx.second() - 1,
                                         
                                         cidx.third());
                    
                    if (genfunc::addValidAndUniqueTriple(recvec, t10idx)) nterms++;
                    
                    // (0 |A(0)| b - 1)^(m+1) term
                    
                    CThreeIndexes t11idx(cidx.first(),  cidx.second() - 1,
                                         
                                         cidx.third() + 1);
                    
                    if (genfunc::addValidAndUniqueTriple(recvec, t11idx)) nterms++;
                    
                    // (0 |A(0)| b - 2)^(m) term
                    
                    CThreeIndexes t20idx(cidx.first(),  cidx.second() - 2,
                                         
                                         cidx.third());
                    
                    if (genfunc::addValidAndUniqueTriple(recvec, t20idx)) nterms++;
                    
                    // (0 |A(0)| b - 2)^(m+1) term
                    
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
CNuclearPotentialIntegralsDriver::_getIndexesForRecursionPattern(      std::vector<int32_t>& recIndexes,
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
