//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef DensityGridDriver_hpp
#define DensityGridDriver_hpp

#include <cstdint>
#include <vector>

#include "mpi.h"

#include "ExecMode.hpp"
#include "GtoContainer.hpp"
#include "MolecularBasis.hpp"
#include "MolecularGrid.hpp"
#include "Molecule.hpp"
#include "SphericalMomentum.hpp"
#include "VecMemBlocks.hpp"
#include "XCFuncType.hpp"
#include "AODensityMatrix.hpp"
#include "DensityGrid.hpp"

/**
 Class CDensityGridDriver generates density grid for usage in numerical
 integration.

 @author Z. Rinkevicius
 */
class CDensityGridDriver
{
    /**
     The rank of associated local MPI process.
     */
    int32_t _locRank;

    /**
     The total number of local MPI processes.
     */
    int32_t _locNodes;

    /**
     The MPI communicator.
     */
    MPI_Comm _locComm;

    /**
     The threshold of density screening.
     */
    double _thresholdOfDensity;

    /**
     The execution mode of grid driver object.
     */
    execmode _runMode;

    /**
     Creates density grid on each MPI node within domain of MPI communicator.
     Density grid points are generated using only CPUs.

     @param densityGrid the density grid object.
     @param aoDensityMatrix the AO density matrices.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param molecularGrid the distributed molecular grid.
     @param xcFunctional the exchange-correlation functional type.
     */
    void _genDensityGridOnCPU(      CDensityGrid&     densityGrid,
                              const CAODensityMatrix& aoDensityMatrix,
                              const CMolecule&        molecule,
                              const CMolecularBasis&  basis,
                              const CMolecularGrid&   molecularGrid,
                              const xcfun             xcFunctional);
    
    /**
     Creates density grid for spin-restricted LDA case.
     
     @param densityGrid the density grid object.
     @param aoDensityMatrix the AO density matrices.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param molecularGrid the distributed molecular grid.
     */
    void _genRestrictedDensityForLDA(      CDensityGrid&     densityGrid,
                                     const CAODensityMatrix& aoDensityMatrix,
                                     const CMolecule&        molecule,
                                     const CMolecularBasis&  basis,
                                     const CMolecularGrid&   molecularGrid) const;
    
    
    /**
     Creates density grid for spin-restricted GGA case.
     
     @param densityGrid the density grid object.
     @param aoDensityMatrix the AO density matrices.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param molecularGrid the distributed molecular grid.
     */
    void _genRestrictedDensityForGGA(      CDensityGrid&     densityGrid,
                                     const CAODensityMatrix& aoDensityMatrix,
                                     const CMolecule&        molecule,
                                     const CMolecularBasis&  basis,
                                     const CMolecularGrid&   molecularGrid) const;
    
    /**
      Generates batch of spin restricted density grid points for LDA case.

     @param densityGrid the pointer to density grid object.
     @param aoDensityMatrix the AO density matrix.
     @param gtoContainer the GTOs container.
     @param gridCoordinatesX the vector of Cartesian X coordinates of grid
     points.
     @param gridCoordinatesY the vector of Cartesian Y coordinates of grid
     points.
     @param gridCoordinatesZ the vector of Cartesian Y coordinates of grid
     points.
     @param gridOffset the batch offset in vector grid points.
     @param nGridPoints the number of grid points in batch.
     */
    void _genBatchOfRestrictedDensityGridPointsForLDA(      CDensityGrid*     densityGrid,
                                                      const CAODensityMatrix* aoDensityMatrix,
                                                      const CGtoContainer*    gtoContainer,
                                                      const double*           gridCoordinatesX,
                                                      const double*           gridCoordinatesY,
                                                      const double*           gridCoordinatesZ,
                                                      const int32_t           gridOffset,
                                                      const int32_t           nGridPoints) const;
    
    /**
     Generates batch of spin restricted density grid points for GGA case.
     
     @param densityGrid the pointer to density grid object.
     @param aoDensityMatrix the AO density matrix.
     @param gtoContainer the GTOs container.
     @param gridCoordinatesX the vector of Cartesian X coordinates of grid
     points.
     @param gridCoordinatesY the vector of Cartesian Y coordinates of grid
     points.
     @param gridCoordinatesZ the vector of Cartesian Y coordinates of grid
     points.
     @param gridOffset the batch offset in vector grid points.
     @param nGridPoints the number of grid points in batch.
     */
    void _genBatchOfRestrictedDensityGridPointsForGGA(      CDensityGrid*     densityGrid,
                                                      const CAODensityMatrix* aoDensityMatrix,
                                                      const CGtoContainer*    gtoContainer,
                                                      const double*           gridCoordinatesX,
                                                      const double*           gridCoordinatesY,
                                                      const double*           gridCoordinatesZ,
                                                      const int32_t           gridOffset,
                                                      const int32_t           nGridPoints) const;
    
    
    /**
     Distributes spin-restriced density values into density grid.

     @param densityGrid the pointer to density grid object.
     @param aoDensityMatrix the AO density matrix.
     @param gtoValues the GTOs values buffer.
     @param gridOffset the offset of grid points batch in molecular grid.
     @param gridBlockPosition the position of grid block in GTOs values grid.
     @param nGridPoints the number of grid points in grid points batch.
     */
    void _distRestrictedDensityValuesForLDA(      CDensityGrid*        densityGrid,
                                            const CAODensityMatrix*    aoDensityMatrix,
                                            const CMemBlock2D<double>& gtoValues,
                                            const int32_t              gridOffset,
                                            const int32_t              gridBlockPosition,
                                            const int32_t              nGridPoints) const;
    
    /**
     Distributes spin-restriced density values into density grid.
     
     @param densityGrid the pointer to density grid object.
     @param aoDensityMatrix the AO density matrix.
     @param gtoValues the GTOs values buffer.
     @param gtoValuesX the GTOs gradient along X axis values buffer.
     @param gtoValuesY the GTOs gradient along Y axis values buffer.
     @param gtoValuesZ the GTOs gradient along Z axis values buffer.
     @param gridOffset the offset of grid points batch in molecular grid.
     @param gridBlockPosition the position of grid block in GTOs values grid.
     @param nGridPoints the number of grid points in grid points batch.
     */
    void _distRestrictedDensityValuesForGGA(      CDensityGrid*        densityGrid,
                                            const CAODensityMatrix*    aoDensityMatrix,
                                            const CMemBlock2D<double>& gtoValues,
                                            const CMemBlock2D<double>& gtoValuesX,
                                            const CMemBlock2D<double>& gtoValuesY,
                                            const CMemBlock2D<double>& gtoValuesZ,
                                            const int32_t              gridOffset,
                                            const int32_t              gridBlockPosition,
                                            const int32_t              nGridPoints) const;
    
    /**
     Gets size of block in grid batch.

     @return the size of block in grid batch.
     */
    int32_t _getSizeOfBlock() const;

   public:
    /**
     Creates a density grid driver object using MPI info.

     @param comm the MPI communicator.
     */
    CDensityGridDriver(MPI_Comm comm);

    /**
     Destroys a grid driver object.
     */
    ~CDensityGridDriver();

    /**
     Generates partitioned density grid for given molecule and type of
     exchange-correlation functional. Density grid generation is distributed
     within domain of MPI communicator.

     @param aoDensityMatrix the AO density matrix.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param molecularGrid the molecular grid.
     @param xcFunctional the type of exchange-correlation functional.
     @return the density grid object.
     */
    CDensityGrid generate(const CAODensityMatrix& aoDensityMatrix,
                          const CMolecule&        molecule,
                          const CMolecularBasis&  basis,
                          const CMolecularGrid&   molecularGrid,
                          const xcfun             xcFunctional);
};

#endif /* DensityGridDriver_hpp */
