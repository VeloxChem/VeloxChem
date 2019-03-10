//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ElectricFieldIntegralsDriver.hpp"

#include "AngularMomentum.hpp"
#include "GenFunc.hpp"
#include "BoysFunction.hpp"
#include "OneIntsFunc.hpp"

CElectricFieldIntegralsDriver::CElectricFieldIntegralsDriver(const int32_t  globRank,
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

CElectricFieldIntegralsDriver::~CElectricFieldIntegralsDriver()
{
}

CElectricFieldMatrix
CElectricFieldIntegralsDriver::compute(const CMolecule&       molecule,
                                       const CMolecularBasis& basis,
                                       const double           coordinateX,
                                       const double           coordinateY,
                                       const double           coordinateZ,
                                             MPI_Comm         comm) const
{
    CElectricFieldMatrix efieldmat;
    
    if (_locRank == mpi::master())
    {
        // set up GTOs container
        
        CGtoContainer bracontr(molecule, basis);
        
        // set up poinr dipoles data
        
        auto dipoles = CMemBlock2D<double>({1.0, 1.0, 1.0}, 1, 3);
        
        auto coords  = CMemBlock2D<double>({coordinateX, coordinateY, coordinateZ},
                                           1, 3);
        
        // compute electric field integrals
        
        efieldmat = _compElectricFieldIntegrals(&dipoles, &coords, &bracontr,
                                                &bracontr);
    }
    
    return efieldmat;
}

CElectricFieldMatrix
CElectricFieldIntegralsDriver::compute(const CMolecule&           molecule,
                                       const CMolecularBasis&     basis,
                                       const CMemBlock2D<double>* dipoles,
                                       const CMemBlock2D<double>* coordinates,
                                             MPI_Comm             comm) const
{
    CElectricFieldMatrix efieldmat;
    
    if (_locRank == mpi::master())
    {
        // set up GTOs container
        
        CGtoContainer bracontr(molecule, basis);
        
        // compute electric field integrals
        
        efieldmat = _compElectricFieldIntegrals(dipoles, coordinates,
                                                &bracontr, &bracontr);
    }
    
    return efieldmat;
}

CElectricFieldMatrix
CElectricFieldIntegralsDriver::compute(const CMolecule&       molecule,
                                       const CMolecularBasis& braBasis,
                                       const CMolecularBasis& ketBasis,
                                       const double           coordinateX,
                                       const double           coordinateY,
                                       const double           coordinateZ,
                                             MPI_Comm         comm) const
{
    CElectricFieldMatrix efieldmat;
    
    if (_locRank == mpi::master())
    {
        // set up GTOs container
        
        CGtoContainer bracontr(molecule, braBasis);
        
        CGtoContainer ketcontr(molecule, ketBasis);
        
        // set up point dipoles data
        
        auto dipoles = CMemBlock2D<double>({1.0, 1.0, 1.0}, 1, 3);
        
        auto coords  = CMemBlock2D<double>({coordinateX, coordinateY, coordinateZ},
                                           1, 3);
        
        // compute electric field integrals
        
        efieldmat = _compElectricFieldIntegrals(&dipoles, &coords, &bracontr,
                                                 &ketcontr);
    }
    
    return efieldmat;
}

CElectricFieldMatrix
CElectricFieldIntegralsDriver::compute(const CMolecule&           molecule,
                                       const CMolecularBasis&     braBasis,
                                       const CMolecularBasis&     ketBasis,
                                       const CMemBlock2D<double>* dipoles,
                                       const CMemBlock2D<double>* coordinates,
                                             MPI_Comm             comm) const
{
    CElectricFieldMatrix efieldmat;
    
    if (_locRank == mpi::master())
    {
        // set up GTOs container
        
        CGtoContainer bracontr(molecule, braBasis);
        
        CGtoContainer ketcontr(molecule, ketBasis);
        
        // compute electric field integrals
        
        efieldmat = _compElectricFieldIntegrals(dipoles, coordinates, &bracontr,
                                                &ketcontr);
    }
    
    return efieldmat;
}

CElectricFieldMatrix
CElectricFieldIntegralsDriver::compute(const CMolecule&       braMolecule,
                                       const CMolecule&       ketMolecule,
                                       const CMolecularBasis& basis,
                                       const double           coordinateX,
                                       const double           coordinateY,
                                       const double           coordinateZ,
                                             MPI_Comm         comm) const
{
  
    CElectricFieldMatrix efieldmat;
    
    if (_locRank == mpi::master())
    {
        // set up GTOs container
        
        CGtoContainer bracontr(braMolecule, basis);
        
        CGtoContainer ketcontr(ketMolecule, basis);
        
        // set up point dipoles data
        
        auto dipoles = CMemBlock2D<double>({1.0, 1.0, 1.0}, 1, 3);
        
        auto coords  = CMemBlock2D<double>({coordinateX, coordinateY, coordinateZ},
                                           1, 3);
        
        // compute electric field integrals
        
        efieldmat = _compElectricFieldIntegrals(&dipoles, &coords, &bracontr,
                                                &ketcontr);
    }
    
    return efieldmat;
}

CElectricFieldMatrix
CElectricFieldIntegralsDriver::compute(const CMolecule&           braMolecule,
                                       const CMolecule&           ketMolecule,
                                       const CMolecularBasis&     basis,
                                       const CMemBlock2D<double>* dipoles,
                                       const CMemBlock2D<double>* coordinates,
                                             MPI_Comm             comm) const
{
    
    CElectricFieldMatrix efieldmat;
    
    if (_locRank == mpi::master())
    {
        // set up GTOs container
        
        CGtoContainer bracontr(braMolecule, basis);
        
        CGtoContainer ketcontr(ketMolecule, basis);
        
        // compute electric field integrals
        
        efieldmat = _compElectricFieldIntegrals(dipoles, coordinates, &bracontr,
                                                &ketcontr);
    }
    
    return efieldmat;
}

CElectricFieldMatrix
CElectricFieldIntegralsDriver::compute(const CMolecule&       braMolecule,
                                       const CMolecule&       ketMolecule,
                                       const CMolecularBasis& braBasis,
                                       const CMolecularBasis& ketBasis,
                                       const double           coordinateX,
                                       const double           coordinateY,
                                       const double           coordinateZ,
                                             MPI_Comm         comm) const
{
    CElectricFieldMatrix efieldmat;
    
    if (_locRank == mpi::master())
    {
        // set up GTOs container
        
        CGtoContainer bracontr(braMolecule, braBasis);
        
        CGtoContainer ketcontr(ketMolecule, ketBasis);
        
        // set up point dipoles data
        
        auto dipoles = CMemBlock2D<double>({1.0, 1.0, 1.0}, 1, 3);
        
        auto coords  = CMemBlock2D<double>({coordinateX, coordinateY, coordinateZ},
                                           1, 3);
        
        // compute nuclear potential integrals
        
        efieldmat = _compElectricFieldIntegrals(&dipoles, &coords, &bracontr,
                                                &ketcontr);
    }
    
    return efieldmat;
}

CElectricFieldMatrix
CElectricFieldIntegralsDriver::compute(const CMolecule&           braMolecule,
                                       const CMolecule&           ketMolecule,
                                       const CMolecularBasis&     braBasis,
                                       const CMolecularBasis&     ketBasis,
                                       const CMemBlock2D<double>* dipoles,
                                       const CMemBlock2D<double>* coordinates,
                                             MPI_Comm             comm) const
{
    CElectricFieldMatrix efieldmat;
    
    if (_locRank == mpi::master())
    {
        // set up GTOs container
        
        CGtoContainer bracontr(braMolecule, braBasis);
        
        CGtoContainer ketcontr(ketMolecule, ketBasis);
        
        // compute nuclear potential integrals
        
        efieldmat = _compElectricFieldIntegrals(dipoles, coordinates, &bracontr,
                                                &ketcontr);
    }
    
    return efieldmat;
}

void
CElectricFieldIntegralsDriver::compute(         double*              intsBatchX,
                                                double*              intsBatchY,
                                                double*              intsBatchZ,
                                          const CMemBlock2D<double>* dipoles,
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
    
    // set up distribution pattern
    
    COneIntsDistribution distx(intsBatchX, nrow, ncol, dist1e::batch);
    
    COneIntsDistribution disty(intsBatchY, nrow, ncol, dist1e::batch);
    
    COneIntsDistribution distz(intsBatchZ, nrow, ncol, dist1e::batch);
    
    // compute electric field integrals
    
    _compElectricFieldForGtoBlocks(&distx, &disty, &distz, dipoles, coordinates,
                                   braGtoBlock, ketGtoBlock);
}

CElectricFieldMatrix
CElectricFieldIntegralsDriver::_compElectricFieldIntegrals(const CMemBlock2D<double>* dipoles,
                                                           const CMemBlock2D<double>* coordinates,
                                                           const CGtoContainer*       braGtoContainer,
                                                           const CGtoContainer*       ketGtoContainer) const
{
    // check if GTOs containers are same on bra and ket sides
    
    auto symbk = ((*braGtoContainer) == (*ketGtoContainer));
    
    // determine dimensions of overlap matrix
    
    auto nrow = braGtoContainer->getNumberOfAtomicOrbitals();
    
    auto ncol = ketGtoContainer->getNumberOfAtomicOrbitals();
    
    // allocate dense matrix for electric field integrals
    
    CDenseMatrix efxmat(nrow, ncol);
    
    CDenseMatrix efymat(nrow, ncol);
    
    CDenseMatrix efzmat(nrow, ncol);
    
    // set up distributio pattern
    
    dist1e dstyp = (symbk) ? dist1e::symsq : dist1e::rect;
    
    COneIntsDistribution* distpatx = new COneIntsDistribution(efxmat.values(),
                                                              nrow, ncol, dstyp);
    
    COneIntsDistribution* distpaty = new COneIntsDistribution(efymat.values(),
                                                              nrow, ncol, dstyp);
    
    COneIntsDistribution* distpatz = new COneIntsDistribution(efzmat.values(),
                                                              nrow, ncol, dstyp);
    
    // compute electric field integral blocks
    
    #pragma omp parallel shared(braGtoContainer, ketGtoContainer, dipoles,\
                                coordinates, distpatx, distpaty, distpatz, symbk)
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
                        
                        _compElectricFieldForGtoBlocks(distpatx, distpaty,
                                                       distpatz, dipoles,
                                                       coordinates, bgtos,
                                                       kgtos);
                    }
                }
            }
        }
    }
    
    // deallocate distribution pattern
    
    delete distpatx;
    
    delete distpaty;
    
    delete distpatz;
    
    return CElectricFieldMatrix(efxmat, efymat, efzmat);
}

void
CElectricFieldIntegralsDriver::_compElectricFieldForGtoBlocks(      COneIntsDistribution* distPatternX,
                                                                    COneIntsDistribution* distPatternY,
                                                                    COneIntsDistribution* distPatternZ,
                                                              const CMemBlock2D<double>*  dipoles,
                                                              const CMemBlock2D<double>*  coordinates,
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
    
    // copy dipoles and their coordinates
    
    auto dipvalues = *dipoles;
    
    auto dipcoords = *coordinates;
    
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
    
    // allocate primitives integrals buffer
    
    CMemBlock2D<double> pbuffer(pdim, nblk);
    
    // allocate accumulation buffer
    
    auto ncart = angmom::to_CartesianComponents(bang, kang);
    
    CMemBlock2D<double> accbuffer(pdim, 3 * ncart * pmax);
    
    // set up contracted GTOs dimensions
    
    auto kdim = ketgtos.getNumberOfContrGtos();
    
    // allocate contracted Cartesian integrals buffer
    
    CMemBlock2D<double> cartbufferx(kdim, ncart);
    
    CMemBlock2D<double> cartbuffery(kdim, ncart);
    
    CMemBlock2D<double> cartbufferz(kdim, ncart);
    
    // allocate contracted spherical integrals buffer
    
    auto nspher = angmom::to_SphericalComponents(bang, kang);
    
    CMemBlock2D<double> spherbufferx(kdim, nspher);
    
    CMemBlock2D<double> spherbuffery(kdim, nspher);
    
    CMemBlock2D<double> spherbufferz(kdim, nspher);
    
    // set up indexes for contraction
    
    auto pidx = genfunc::findQuadrupleIndex(recidx, recvec, {bang, kang, 0, 0});
    
    auto spos = braGtoBlock.getStartPositions();
    
    auto epos = braGtoBlock.getEndPositions();
    
    // initialize Boys function evaluator
    
    auto bord = genfunc::maxOrderOfPair(recvec, 0, 0);
    
    CBoysFunction bftab(bord);
    
    CMemBlock<double> bargs(pdim);
    
    CMemBlock2D<double> bvals(pdim, bord + 1);
    
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
        
        for (int32_t j = 0; j < dipvalues.size(0); j++)
        {
            // compute distances: R(PC) = P - C
            
            intsfunc::compDistancesPC(rpc, rp, dipcoords, bragtos, ketgtos, i, j);
            
            // compute primitive integrals
            
            //_compPrimNuclearPotentialInts(pbuffer, recvec, recidx, bftab, bargs,
            //                              bvals, bord, rfacts, rab, rpa, rpb,
            //                              rpc, bragtos, ketgtos, i);
            
            // add scaled contribution to accumulation buffer
            
            _addPointDipoleContribution(accbuffer, pbuffer, pidx, dipvalues,
                                        bragtos, ketgtos, i, j);
        }
        
        // contract primitive overlap integrals
        
        genfunc::contract(cartbufferx, accbuffer, 0, bragtos, ketgtos, i);
        
        genfunc::contract(cartbuffery, accbuffer, poff, bragtos, ketgtos, i);
        
        genfunc::contract(cartbufferz, accbuffer, 2 * poff, bragtos, ketgtos, i);
        
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

CVecFourIndexes
CElectricFieldIntegralsDriver::_getRecursionPattern(const CGtoBlock& braGtoBlock,
                                                    const CGtoBlock& ketGtoBlock) const
{
    // set up angular momentum
    
    auto bang = braGtoBlock.getAngularMomentum();
    
    auto kang = ketGtoBlock.getAngularMomentum();
    
    // set up recursion buffer
    
    CVecFourIndexes recvec;
    
    recvec.reserve((bang + 1) * (kang + 1));
    
    // set up indexing counters
    
    int32_t spos = 0;
    
    int32_t epos = 1;
    
    // set up initial state of recursion buffer
    
    recvec.push_back(CFourIndexes(bang, kang, 0, 0));
    
    while (true)
    {
        // internal new recursion terms counter
        
        int32_t nterms = 0;
        
        // generate bra and ket Obara-Saika recursion terms
        
        for (int32_t i = spos; i < epos; i++)
        {
            CFourIndexes cidx(recvec[i]);
            
            if (cidx.fourth() == 0)
            {
                // electric field recursion
                
                if (cidx.first() != 0)
                {
                    // general recursion for bra and ket sides
                
                    // (a - 1 |A(1)| b)^(m) term
                
                    CFourIndexes t10idx(cidx.first() - 1,  cidx.second(),
                                        
                                        cidx.third(), 0);
                
                    if (genfunc::addValidAndUniqueQuadruple(recvec, t10idx)) nterms++;
                
                    // (a - 1 |A(1)| b)^(m+1) term
                
                    CFourIndexes t11idx(cidx.first() - 1,  cidx.second(),
                                        
                                        cidx.third() + 1, 0);
                
                    if (genfunc::addValidAndUniqueQuadruple(recvec, t11idx)) nterms++;
                
                    // (a - 2 |A(1)| b)^(m) term
                
                    CFourIndexes t20idx(cidx.first() - 2,  cidx.second(),
                                        
                                        cidx.third(), 0);
                
                    if (genfunc::addValidAndUniqueQuadruple(recvec, t20idx)) nterms++;
                
                    // (a - 2 |A(1)| b)^(m+1) term
                
                    CFourIndexes t21idx(cidx.first() - 2,  cidx.second(),
                                        
                                        cidx.third() + 1, 0);
                
                    if (genfunc::addValidAndUniqueQuadruple(recvec, t21idx)) nterms++;
                
                    // (a - 1 |A(1)| b - 1)^(m) term
                
                    CFourIndexes tk0idx(cidx.first() - 1,  cidx.second() - 1,
                                        
                                        cidx.third(), 0);
                
                    if (genfunc::addValidAndUniqueQuadruple(recvec, tk0idx)) nterms++;
                
                    // (a - 1 |A(1)| b - 1)^(m+1) term
                
                    CFourIndexes tk1idx(cidx.first() - 1,  cidx.second() - 1,
                                        
                                        cidx.third() + 1, 0);
                
                    if (genfunc::addValidAndUniqueQuadruple(recvec, tk1idx)) nterms++;
                
                    // (a - 1 |A(0)| b)^(m+1) term
                
                    CFourIndexes s10idx(cidx.first() - 1,  cidx.second(),
                                        
                                        cidx.third() + 1, 1);
                
                    if (genfunc::addValidAndUniqueQuadruple(recvec, s10idx)) nterms++;
                }
                else
                {
                    // simplified recursion for ket sides
                    
                    // (0 |A(1)| b - 1)^(m) term
                    
                    CFourIndexes t10idx(cidx.first(),  cidx.second() - 1,
                                        
                                        cidx.third(), 0);
                    
                    if (genfunc::addValidAndUniqueQuadruple(recvec, t10idx)) nterms++;
                    
                    // (0 |A(1)| b - 1)^(m+1) term
                    
                    CFourIndexes t11idx(cidx.first(),  cidx.second() - 1,
                                        
                                        cidx.third() + 1, 0);
                    
                    if (genfunc::addValidAndUniqueQuadruple(recvec, t11idx)) nterms++;
                    
                    // (0 |A(1)| b - 2)^(m) term
                    
                    CFourIndexes t20idx(cidx.first(),  cidx.second() - 2,
                                         
                                         cidx.third(), 0);
                    
                    if (genfunc::addValidAndUniqueQuadruple(recvec, t20idx)) nterms++;
                    
                    // (0 |A(1)| b - 2)^(m+1) term
                    
                    CFourIndexes t21idx(cidx.first(),  cidx.second() - 2,
                                         
                                        cidx.third() + 1, 0);
                    
                    if (genfunc::addValidAndUniqueQuadruple(recvec, t21idx)) nterms++;
                    
                    // (0 |A(0)| b - 1)^(m+1) term
                    
                    CFourIndexes s10idx(cidx.first(),  cidx.second() - 1,
                                        
                                        cidx.third() + 1 , 1);
                    
                    if (genfunc::addValidAndUniqueQuadruple(recvec, s10idx)) nterms++;
                }
            }
            else
            {
                // nuclear repulsion recursion
                
                if (cidx.first() != 0)
                {
                    // general recursion for bra and ket sides
                    
                    // (a - 1 |A(0)| b)^(m) term
                    
                    CFourIndexes s10idx(cidx.first() - 1,  cidx.second(),
                                        
                                        cidx.third(), 1);
                    
                    if (genfunc::addValidAndUniqueQuadruple(recvec, s10idx)) nterms++;
                    
                    // (a - 1 |A(0)| b)^(m+1) term
                    
                    CFourIndexes s11idx(cidx.first() - 1,  cidx.second(),
                                        
                                        cidx.third() + 1, 1);
                    
                    if (genfunc::addValidAndUniqueQuadruple(recvec, s11idx)) nterms++;
                    
                    // (a - 2 |A(0)| b)^(m) term
                    
                    CFourIndexes s20idx(cidx.first() - 2,  cidx.second(),
                                        
                                        cidx.third(), 1);
                    
                    if (genfunc::addValidAndUniqueQuadruple(recvec, s20idx)) nterms++;
                    
                    // (a - 2 |A(0)| b)^(m+1) term
                    
                    CFourIndexes s21idx(cidx.first() - 2,  cidx.second(),
                                        
                                        cidx.third() + 1, 1);
                    
                    if (genfunc::addValidAndUniqueQuadruple(recvec, s21idx)) nterms++;
                    
                    // (a - 1 |A(0)| b - 1)^(m) term
                    
                    CFourIndexes sk0idx(cidx.first() - 1,  cidx.second() - 1,
                                        
                                        cidx.third(), 1);
                    
                    if (genfunc::addValidAndUniqueQuadruple(recvec, sk0idx)) nterms++;
                    
                    // (a - 1 |A(0)| b - 1)^(m+1) term
                    
                    CFourIndexes sk1idx(cidx.first() - 1,  cidx.second() - 1,
                                        
                                        cidx.third() + 1, 1);
                    
                    if (genfunc::addValidAndUniqueQuadruple(recvec, sk1idx)) nterms++;
                }
                else
                {
                    // simplified recursion for ket sides
                    
                    // (0 |A(0)| b - 1)^(m) term
                    
                    CFourIndexes s10idx(cidx.first(),  cidx.second() - 1,
                                        
                                        cidx.third(), 1);
                    
                    if (genfunc::addValidAndUniqueQuadruple(recvec, s10idx)) nterms++;
                    
                    // (0 |A(0)| b - 1)^(m+1) term
                    
                    CFourIndexes s11idx(cidx.first(),  cidx.second() - 1,
                                        
                                        cidx.third() + 1, 1);
                    
                    if (genfunc::addValidAndUniqueQuadruple(recvec, s11idx)) nterms++;
                    
                    // (0 |A(0)| b - 2)^(m) term
                    
                    CFourIndexes s20idx(cidx.first(),  cidx.second() - 2,
                                        
                                        cidx.third(), 1);
                    
                    if (genfunc::addValidAndUniqueQuadruple(recvec, s20idx)) nterms++;
                    
                    // (0 |A(0)| b - 2)^(m+1) term
                    
                    CFourIndexes s21idx(cidx.first(),  cidx.second() - 2,
                                        
                                        cidx.third() + 1, 1);
                    
                    if (genfunc::addValidAndUniqueQuadruple(recvec, s21idx)) nterms++;
                }
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
CElectricFieldIntegralsDriver::_getIndexesForRecursionPattern(      std::vector<int32_t>& recIndexes,
                                                              const CVecFourIndexes&     recPattern,
                                                              const int32_t              maxPrimGtos) const
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
        
        if (recPattern[i].fourth() == 0)
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

void
CElectricFieldIntegralsDriver::_addPointDipoleContribution(      CMemBlock2D<double>& accBuffer,
                                                           const CMemBlock2D<double>& primBuffer,
                                                           const int32_t              primIndex,
                                                           const CMemBlock2D<double>& dipoles,
                                                           const CGtoBlock&           braGtoBlock,
                                                           const CGtoBlock&           ketGtoBlock,
                                                           const int32_t              iContrGto,
                                                           const int32_t              iPointDipole) const
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
    
    for (int32_t i = 0; i < 3; i++)
    {
        // set up point dipole factor
        
        auto fact = (dipoles.data(i))[iPointDipole];
        
        // set up electric field component offset
        
        auto poff = i * bdim * ncart;
        
        for (int32_t j = 0; j < bdim; j++)
        {
            for (int32_t k = 0; k < ncart; k++)
            {
                auto abuf = accBuffer.data(poff + j * ncart + k);
                
                auto pbuf = primBuffer.data(primIndex + poff + j * ncart + k);
                
                #pragma omp simd aligned(abuf, pbuf: VLX_ALIGN)
                for (int32_t l = 0; l < nprim; l++)
                {
                    abuf[l] += fact * pbuf[l];
                }
            }
        }
    }
}
