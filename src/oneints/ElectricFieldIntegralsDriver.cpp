//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ElectricFieldIntegralsDriver.hpp"

#include "AngularMomentum.hpp"

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
                // electric field integrals
                
                
            }
            else
            {
                // nuclear repulsion integrals
                
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
