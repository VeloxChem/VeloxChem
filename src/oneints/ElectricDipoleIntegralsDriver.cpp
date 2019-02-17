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
#include "ElectricDipoleRecFunc.hpp"

CElectricDipoleIntegralsDriver::CElectricDipoleIntegralsDriver(const int32_t  globRank,
                                                               const int32_t  globNodes,
                                                                     MPI_Comm comm)

    : _globRank(globRank)

    , _globNodes(globNodes)

    , _isLocalMode(false)

    , _xOrigin(0.0)

    , _yOrigin(0.0)

    , _zOrigin(0.0)
{
    _locRank  = mpi::rank(comm);
    
    _locNodes = mpi::nodes(comm);
    
    _isLocalMode = !mpi::compare(comm, MPI_COMM_WORLD);
}

CElectricDipoleIntegralsDriver::~CElectricDipoleIntegralsDriver()
{
    
}

void
CElectricDipoleIntegralsDriver::setElectricDipoleOrigin(const double xOrigin,
                                                        const double yOrigin,
                                                        const double zOrigin)
{
    _xOrigin = xOrigin;
    
    _yOrigin = yOrigin;
    
    _zOrigin = zOrigin; 
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
    
    _compElectricDipoleForGtoBlocks(&distx, &disty, &distz, _xOrigin, _yOrigin,
                                    _zOrigin, braGtoBlock, ketGtoBlock);
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
    
    COneIntsDistribution* distpaty = new COneIntsDistribution(dipymat.values(),
                                                              nrow, ncol, dstyp);
    
    COneIntsDistribution* distpatz = new COneIntsDistribution(dipzmat.values(),
                                                              nrow, ncol, dstyp);
    
    // copy origin coordinates
    
    auto origx = _xOrigin;
    
    auto origy = _yOrigin;
    
    auto origz = _zOrigin;
    
    // compute overlap integral blocks
    
    #pragma omp parallel shared(braGtoContainer, ketGtoContainer, distpatx,\
                                distpaty, distpatz, origx, origy, origz, symbk)
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
                                                        distpatz, origx, origy,
                                                        origz, bgtos, kgtos);
                    }
                }
            }
        }
    }
    
    // deallocate distribution pattern
    
    delete distpatx;
    
    delete distpaty;
    
    delete distpatz;
    
    return CElectricDipoleMatrix(dipxmat, dipymat, dipzmat, _xOrigin, _yOrigin,
                                 _zOrigin);
}

void
CElectricDipoleIntegralsDriver::_compElectricDipoleForGtoBlocks(      COneIntsDistribution* distPatternX,
                                                                      COneIntsDistribution* distPatternY,
                                                                      COneIntsDistribution* distPatternZ,
                                                                const double                xOrigin,
                                                                const double                yOrigin,
                                                                const double                zOrigin, 
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
    
    // allocate primitives integrals buffer
    
    CMemBlock2D<double> pbuffer(pdim, nblk);
    
    // set up contracted GTOs dimensions
    
    auto kdim = ketgtos.getNumberOfContrGtos();
    
    // allocate contracted Cartesian integrals buffer
    
    auto ncart = angmom::to_CartesianComponents(bang, kang);
    
    CMemBlock2D<double> cartbufferx(kdim, ncart);
    
    CMemBlock2D<double> cartbuffery(kdim, ncart);
    
    CMemBlock2D<double> cartbufferz(kdim, ncart);
    
    // allocate contracted spherical integrals buffer
    
    auto nspher = angmom::to_SphericalComponents(bang, kang);
    
    CMemBlock2D<double> spherbufferx(kdim, nspher);
    
    CMemBlock2D<double> spherbuffery(kdim, nspher);
    
    CMemBlock2D<double> spherbufferz(kdim, nspher);
    
    // set up indexes for contraction
    
    auto pidx = genfunc::findTripleIndex(recidx, recvec, {bang, kang, 0});
    
    auto spos = braGtoBlock.getStartPositions();
    
    auto epos = braGtoBlock.getEndPositions();
    
    // determine bra and ket sides symmetry
    
    bool symbk = (bragtos == ketgtos);

    for (int32_t i = 0; i < bragtos.getNumberOfContrGtos(); i++)
    {
        // compute bra dimensions and shift
        
        auto bdim = epos[i] - spos[i];
        
        auto poff = bdim * ncart;
        
        // compute distances: R(AB) = A - B
        
        intsfunc::compDistancesAB(rab, bragtos, ketgtos, i);
        
        // compute Obara-Saika recursion factors
        
        intsfunc::compFactorsForOverlap(rfacts, bragtos, ketgtos, i);
        
        // compute coordinates of center P
        
        intsfunc::compCoordinatesForP(rp, rfacts, 2, bragtos, ketgtos, i);
        
        // compute distances: R(PA) = P - A
        
        intsfunc::compDistancesPA(rpa, rp, bragtos, ketgtos, i);
        
        // compute distances: R(PB) = P - B
        
        intsfunc::compDistancesPB(rpb, rp, bragtos, ketgtos, i);
        
        // compute distances: R(PC) = P - C
        
        intsfunc::compDistancesPC(rpc, rp, _xOrigin, _yOrigin, _zOrigin,
                                  bragtos, ketgtos, i);
        
        // compute primitive electric dipole integrals
        
        _compPrimElectricDipoleInts(pbuffer, recvec, recidx, rfacts, rab, rpa,
                                    rpb, rpc, bragtos, ketgtos, i);
        
        // contract primitive overlap integrals
        
        genfunc::contract(cartbufferx, pbuffer, pidx, bragtos, ketgtos, i);
        
        genfunc::contract(cartbuffery, pbuffer, pidx + poff, bragtos,
                          ketgtos, i);
        
        genfunc::contract(cartbufferz, pbuffer, pidx + 2 * poff, bragtos,
                          ketgtos, i);
        
        // transform Cartesian to spherical integrals
        
        genfunc::transform(spherbufferx, cartbufferx, bmom, kmom, 0, 0, kdim);
        
        genfunc::transform(spherbuffery, cartbuffery, bmom, kmom, 0, 0, kdim);
        
        genfunc::transform(spherbufferz, cartbufferz, bmom, kmom, 0, 0, kdim);
        
        // add batch of integrals to integrals matrix
        
        distpatx.distribute(spherbufferx, bragtos, ketgtos, symbk, i);
        
        distpaty.distribute(spherbuffery, bragtos, ketgtos, symbk, i);
        
        distpatz.distribute(spherbufferz, bragtos, ketgtos, symbk, i);
    }
}

void
CElectricDipoleIntegralsDriver::_compPrimElectricDipoleInts(      CMemBlock2D<double>&  primBuffer,
                                                            const CVecThreeIndexes&     recPattern,
                                                            const std::vector<int32_t>& recIndexes,
                                                            const CMemBlock2D<double>&  osFactors,
                                                            const CMemBlock2D<double>&  abDistances,
                                                            const CMemBlock2D<double>&  paDistances,
                                                            const CMemBlock2D<double>&  pbDistances,
                                                            const CMemBlock2D<double>&  pcDistances,
                                                            const CGtoBlock&            braGtoBlock,
                                                            const CGtoBlock&            ketGtoBlock,
                                                            const int32_t               iContrGto) const
{
    // compute (s|M|s) integrals
    
    ediprecfunc::compElectricDipoleForSS(primBuffer, recPattern, recIndexes,
                                         osFactors, abDistances, pcDistances,
                                         braGtoBlock, ketGtoBlock, iContrGto);
    
    // compute (s|M|p) integrals
    
    ediprecfunc::compElectricDipoleForSP(primBuffer, recPattern, recIndexes,
                                         osFactors, pbDistances, braGtoBlock,
                                         ketGtoBlock, iContrGto);
    
    // compute (p|M|s) integrals
    
    ediprecfunc::compElectricDipoleForPS(primBuffer, recPattern, recIndexes,
                                         osFactors, paDistances, braGtoBlock,
                                         ketGtoBlock, iContrGto);
    
    // compute (p|M|p) integrals
    
    ediprecfunc::compElectricDipoleForPP(primBuffer, recPattern, recIndexes,
                                         osFactors, paDistances, pbDistances,
                                         braGtoBlock, ketGtoBlock, iContrGto);
    
    // compute (s|M|d) integrals
    
    ediprecfunc::compElectricDipoleForSD(primBuffer, recPattern, recIndexes,
                                         osFactors, pbDistances, braGtoBlock,
                                         ketGtoBlock, iContrGto);
    
    // compute (d|M|s) integrals
    
    ediprecfunc::compElectricDipoleForDS(primBuffer, recPattern, recIndexes,
                                         osFactors, paDistances, braGtoBlock,
                                         ketGtoBlock, iContrGto);
    
    // compute (p|M|d) integrals
    
    ediprecfunc::compElectricDipoleForPD(primBuffer, recPattern, recIndexes,
                                         osFactors, paDistances, pbDistances,
                                         braGtoBlock, ketGtoBlock, iContrGto);
    
    // compute (d|M|p) integrals
    
    ediprecfunc::compElectricDipoleForDP(primBuffer, recPattern, recIndexes,
                                         osFactors, paDistances, braGtoBlock,
                                         ketGtoBlock, iContrGto);
    
    // compute (s|M|f) integrals
    
    ediprecfunc::compElectricDipoleForSF(primBuffer, recPattern, recIndexes,
                                         osFactors, pbDistances, braGtoBlock,
                                         ketGtoBlock, iContrGto);
    
    // compute (f|M|s) integrals
    
    ediprecfunc::compElectricDipoleForFS(primBuffer, recPattern, recIndexes,
                                         osFactors, paDistances, braGtoBlock,
                                         ketGtoBlock, iContrGto);
    
    // compute (d|M|d) integrals
    
    ediprecfunc::compElectricDipoleForDD(primBuffer, recPattern, recIndexes,
                                         osFactors, paDistances, braGtoBlock,
                                         ketGtoBlock, iContrGto);
    
    // compute (p|M|f) integrals
    
    ediprecfunc::compElectricDipoleForPF(primBuffer, recPattern, recIndexes,
                                         osFactors, paDistances, pbDistances,
                                         braGtoBlock, ketGtoBlock, iContrGto);
    
    // compute (f|M|p) integrals
    
    ediprecfunc::compElectricDipoleForFP(primBuffer, recPattern, recIndexes,
                                         osFactors, paDistances, braGtoBlock,
                                         ketGtoBlock, iContrGto);
    
    // compute (s|M|g) integrals
    
    ediprecfunc::compElectricDipoleForSG(primBuffer, recPattern, recIndexes,
                                         osFactors, pbDistances, braGtoBlock,
                                         ketGtoBlock, iContrGto);
    
    // compute (g|M|s) integrals
    
    ediprecfunc::compElectricDipoleForGS(primBuffer, recPattern, recIndexes,
                                         osFactors, paDistances, braGtoBlock,
                                         ketGtoBlock, iContrGto);
    
    // compute (d|M|f) integrals
    
    ediprecfunc::compElectricDipoleForDF(primBuffer, recPattern, recIndexes,
                                         osFactors, paDistances, braGtoBlock,
                                         ketGtoBlock, iContrGto);
    
    // compute (f|M|d) integrals
    
    ediprecfunc::compElectricDipoleForFD(primBuffer, recPattern, recIndexes,
                                         osFactors, paDistances, braGtoBlock,
                                         ketGtoBlock, iContrGto);
    
    // compute (f|M|f) integrals
    
    ediprecfunc::compElectricDipoleForFF(primBuffer, recPattern, recIndexes,
                                         osFactors, paDistances, braGtoBlock,
                                         ketGtoBlock, iContrGto);
    
    // compute (p|M|g) integrals
    
    ediprecfunc::compElectricDipoleForPG(primBuffer, recPattern, recIndexes,
                                         osFactors, paDistances, pbDistances,
                                         braGtoBlock, ketGtoBlock, iContrGto);
    
    // compute (g|M|p) integrals
    
    ediprecfunc::compElectricDipoleForGP(primBuffer, recPattern, recIndexes,
                                         osFactors, paDistances, braGtoBlock,
                                         ketGtoBlock, iContrGto);
    
    // compute (d|M|g) integrals
    
    ediprecfunc::compElectricDipoleForDG(primBuffer, recPattern, recIndexes,
                                         osFactors, paDistances, braGtoBlock,
                                         ketGtoBlock, iContrGto);
    
    // compute (g|M|d) integrals
    
    ediprecfunc::compElectricDipoleForGD(primBuffer, recPattern, recIndexes,
                                         osFactors, paDistances, braGtoBlock,
                                         ketGtoBlock, iContrGto);
    
    // compute (f|M|g) integrals
    
    ediprecfunc::compElectricDipoleForFG(primBuffer, recPattern, recIndexes,
                                         osFactors, paDistances, braGtoBlock,
                                         ketGtoBlock, iContrGto);
    
    // compute (g|M|f) integrals
    
    ediprecfunc::compElectricDipoleForGF(primBuffer, recPattern, recIndexes,
                                         osFactors, paDistances, braGtoBlock,
                                         ketGtoBlock, iContrGto);
    
    // compute (g|M|g) integrals
    
    ediprecfunc::compElectricDipoleForGG(primBuffer, recPattern, recIndexes,
                                         osFactors, paDistances, braGtoBlock,
                                         ketGtoBlock, iContrGto);
    
    // NOTE: add l > 4 recursion here
}

CVecThreeIndexes
CElectricDipoleIntegralsDriver::_getRecursionPattern(const CGtoBlock& braGtoBlock,
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
                // electric dipole recursion
                
                if (cidx.first() != 0)
                {
                    // general recursion for bra and ket sides
                    
                    // (a - 1 |M| b) term
                    
                    CThreeIndexes t1idx(cidx.first() - 1,  cidx.second(), 0);
                    
                    if (genfunc::addValidAndUniqueTriple(recvec, t1idx)) nterms++;
                    
                    // (a - 2 |M| b) term
                    
                    CThreeIndexes t2idx(cidx.first() - 2,  cidx.second(), 0);
                    
                    if (genfunc::addValidAndUniqueTriple(recvec, t2idx)) nterms++;
                    
                    // (a - 1 |M| b - 1) term
                    
                    CThreeIndexes tkidx(cidx.first() - 1,  cidx.second() - 1, 0);
                    
                    if (genfunc::addValidAndUniqueTriple(recvec, tkidx)) nterms++;
                    
                    // (a - 1| b) term
                    
                    CThreeIndexes s0idx(cidx.first() - 1,  cidx.second(), 1);
                    
                    if (genfunc::addValidAndUniqueTriple(recvec, s0idx)) nterms++;
                }
                else
                {
                    // reduced recursion for ket side
                    
                    // (0 |M| b - 1) term
                    
                    CThreeIndexes t1idx(cidx.first(),  cidx.second() - 1, 0);
                    
                    if (genfunc::addValidAndUniqueTriple(recvec, t1idx)) nterms++;
                    
                    // (0 |M| b - 2) term
                    
                    CThreeIndexes t2idx(cidx.first(),  cidx.second() - 2, 0);
                    
                    if (genfunc::addValidAndUniqueTriple(recvec, t2idx)) nterms++;
                    
                    // (0 | b - 1) term
                    
                    CThreeIndexes s0idx(cidx.first(),  cidx.second() - 1, 1);
                    
                    if (genfunc::addValidAndUniqueTriple(recvec, s0idx)) nterms++;
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
   
    // allways compute overlap (s|s)
    
    genfunc::addValidAndUniqueTriple(recvec, CThreeIndexes(0, 0, 1));
    
    return recvec;
}

int32_t
CElectricDipoleIntegralsDriver::_getIndexesForRecursionPattern(      std::vector<int32_t>& recIndexes,
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
        
        auto cblk = maxPrimGtos * angmom::to_CartesianComponents(recPattern[i].first(),
                                                                 recPattern[i].second());
        
        if (recPattern[i].third() == 0)
        {
            nblk += 3 * cblk;
        }
        else
        {
            nblk += cblk;
        }
    }
    
    return nblk;
}

