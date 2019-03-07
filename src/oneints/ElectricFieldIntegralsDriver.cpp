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
    CElectricFieldMatrix efieldmat;
    
    // FIX ME:
    
    return efieldmat;
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
    // FIX ME:
}
