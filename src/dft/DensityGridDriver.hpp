//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

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

/**
 Class CDensityGridDriver generates density grid for usage in numerical
 integration.

 @author Z. Rinkevicius
 */
class CDensityGridDriver
{
    /**
     The rank of associated global MPI process.
     */
    int32_t _globRank;

    /**
     The total number of global MPI processes.
     */
    int32_t _globNodes;

    /**
     The rank of associated local MPI process.
     */
    int32_t _locRank;

    /**
     The total number of local MPI processes.
     */
    int32_t _locNodes;

    /**
     The flag for local execution mode.
     */
    bool _isLocalMode;

    /**
     The threshold of density screening.
     */
    double _thresholdOfDensity;

    /**
     The threshold of primitive GTOs screening.
     */
    double _thresholdOfPrimGTOs;

    /**
     The execution mode of grid driver object.
     */
    execmode _runMode;

    /**
     Creates density grid on each MPI node within domain of MPI communicator.
     Density grid points are generated using only CPUs.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param molGrid the distributed molecular grid.
     @param xcFunctional the exchange-correlation functional type.
     @param comm the MPI communicator.
     */
    void _genDensityGridOnCPU(const CMolecule&       molecule,
                              const CMolecularBasis& basis,
                              const CMolecularGrid&  molGrid,
                              const xcfun            xcFunctional,
                              MPI_Comm               comm);

    /**
     Generates batch of density grid points.

     @param gtoContainer the GTOs container.
     @param gridCoordinatesX the vector of Cartesian X coordinates of grid
            points.
     @param gridCoordinatesY the vector of Cartesian Y coordinates of grid
            points.
     @param gridCoordinatesZ the vector of Cartesian Y coordinates of grid
            points.
     @param gridOffset the batch offset in vector grid points.
     @param nGridPoints the number of grid points in batch.
     @param xcFunctional the exchange-correlation functional type.
     */
    void _genBatchOfDensityGridPoints(const CGtoContainer* gtoContainer,
                                      const double*        gridCoordinatesX,
                                      const double*        gridCoordinatesY,
                                      const double*        gridCoordinatesZ,
                                      const int32_t        gridOffset,
                                      const int32_t        nGridPoints,
                                      const xcfun          xcFunctional);

    /**
     Computes screening factors  vector for primitive Gaussian functions for
     specific grid point using GTOs container.

     @param screenFactors the screening factors vector for each GTOs block
     object in GTOs container.
     @param gtoContainer the GTOs container.
     @param gridCoordinateX the Cartesian X coordinate of grid point.
     @param gridCoordinateY the Cartesian Y coordinate of grid point.
     @param gridCoordinateZ the Cartesian Z coordinate of grid point.
     */
    void _compScreeningFactors(CVecMemBlock<double>& screenFactors,
                               const CGtoContainer&  gtoContainer,
                               const double          gridCoordinateX,
                               const double          gridCoordinateY,
                               const double          gridCoordinateZ);

    /**
     Gets scaling factor of screening factor for specific angular momemntum.

     @param angularMomentum the angular momentum.
     @param maxRadius the maximum radius.
     @return the scaling factor.
     */
    double _getScaleFactor(const int32_t angularMomentum, const double maxRadius) const;

    /**
     Computes vector of distances between grid point and centers of primitive
     Gaussian functions for specific GTOs block.

     @param distances the vector of distances.
     @param gtoContainer the GTOs container.
     @param redDimensions the vector of reduced dimensions.
     @param iGtoBlock the index of GTOs block.
     @param gridCoordinateX the Cartesian X coordinate of grid point.
     @param gridCoordinateY the Cartesian Y coordinate of grid point.
     @param gridCoordinateZ the Cartesian Z coordinate of grid point.
     */
    void _compDistances(CMemBlock2D<double>&        distances,
                        const CGtoContainer&        gtoContainer,
                        const CMemBlock2D<int32_t>& redDimensions,
                        const int32_t               iGtoBlock,
                        const double                gridCoordinateX,
                        const double                gridCoordinateY,
                        const double                gridCoordinateZ) const;

    /**
     Gets number of components required for compupation of all relevant density
     contributions for specific type of exchange-correlation functions.

     @param xcFunctional the exchange-correlations functional type.
     @return the number of components.
     */
    int32_t _getNumberOfXCComponents(const xcfun xcFunctional) const;

    /**
     Computes primitive Gaussian function values for specific GTOs block.

     @param gtoValues the vector of primitive Gaussian functions.
     @param distances the vector of distances.
     @param gtoContainer the GTOs container.
     @param redDimensions the vector of reduced dimensions.
     @param iGtoBlock the index of GTOs block.
     @param xcFunctional the exchange-correlation functional type.
     */
    void _compPrimGtoValues(CMemBlock2D<double>&        gtoValues,
                            const CMemBlock2D<double>&  distances,
                            const CGtoContainer&        gtoContainer,
                            const CMemBlock2D<int32_t>& redDimensions,
                            const int32_t               iGtoBlock,
                            const xcfun                 xcFunctional) const;

    /**
     Computes primitive Gaussian function values for specific GTOs block at
     local density approximation level.

     @param gtoValues the vector of primitive Gaussian functions.
     @param distances the vector of distances.
     @param gtoContainer the GTOs container.
     @param redDimensions the vector of reduced dimensions.
     @param iGtoBlock the index of GTOs block.
     */
    void _compPrimGtoValuesForLDA(CMemBlock2D<double>&        gtoValues,
                                  const CMemBlock2D<double>&  distances,
                                  const CGtoContainer&        gtoContainer,
                                  const CMemBlock2D<int32_t>& redDimensions,
                                  const int32_t               iGtoBlock) const;

    /**
     Computes primitive Gaussian function values for specific GTOs block at
     generalized gradient approximation level.

     @param gtoValues the vector of primitive Gaussian functions.
     @param distances the vector of distances.
     @param gtoContainer the GTOs container.
     @param redDimensions the vector of reduced dimensions.
     @param iGtoBlock the index of GTOs block.
     */
    void _compPrimGtoValuesForGGA(CMemBlock2D<double>&        gtoValues,
                                  const CMemBlock2D<double>&  distances,
                                  const CGtoContainer&        gtoContainer,
                                  const CMemBlock2D<int32_t>& redDimensions,
                                  const int32_t               iGtoBlock) const;

    /**
     Computes primitive Gaussian function values for specific GTOs block at
     meta generalized gradient approximation level.

     @param gtoValues the vector of primitive Gaussian functions.
     @param distances the vector of distances.
     @param gtoContainer the GTOs container.
     @param redDimensions the vector of reduced dimensions.
     @param iGtoBlock the index of GTOs block.
     */
    void _compPrimGtoValuesForMGGA(CMemBlock2D<double>&        gtoValues,
                                   const CMemBlock2D<double>&  distances,
                                   const CGtoContainer&        gtoContainer,
                                   const CMemBlock2D<int32_t>& redDimensions,
                                   const int32_t               iGtoBlock) const;

    /**
     Contracts primitive Gaussian function values for specific GTOs block.

     @param cartGtoValues the vector of contracted Cartesian Gaussian function
            values.
     @param primGtoValues the vector of primitive Gaussian function values.
     @param gtoContainer the GTOs container.
     @param redDimensions the vector of reduced dimensions.
     @param iGtoBlock the index of GTOs block.
     */
    void _contrPrimGtoValues(CMemBlock2D<double>&        cartGtoValues,
                             const CMemBlock2D<double>&  primGtoValues,
                             const CGtoContainer&        gtoContainer,
                             const CMemBlock2D<int32_t>& redDimensions,
                             const int32_t               iGtoBlock) const;

    /**
     Transforms contracted Cartesian GTOs values to spherical GTOs values.

     @param spherGtoValues the spherical GTOs values.
     @param cartGtoValues the Cartesian GTOs values.
     @param spherMomentum the sphericla momentum object.
     @param redDimensions the vector of reduced dimensions.
     @param iGtoBlock the index of GTOs block.
     @param xcFunctional the exchange-correlation functional type.
     */
    void _transContrGtoValues(CMemBlock2D<double>&        spherGtoValues,
                              const CMemBlock2D<double>&  cartGtoValues,
                              const CSphericalMomentum&   spherMomentum,
                              const CMemBlock2D<int32_t>& redDimensions,
                              const int32_t               iGtoBlock,
                              const xcfun                 xcFunctional) const;

   public:
    /**
     Creates a density grid driver object using MPI info.

     @param globRank the the rank of MPI process.
     @param globNodes the total number of MPI processes.
     @param runMode the execution mode.
     @param comm the MPI communicator.
     */
    CDensityGridDriver(const int32_t globRank, const int32_t globNodes, const execmode runMode, MPI_Comm comm);

    /**
     Destroys a grid driver object.
     */
    ~CDensityGridDriver();

    /**
     Generates partitioned density grid for given molecule and type of
     exchange-correlation functional. Density grid generation is distributed
     within domain of MPI communicator.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param molGrid the molecular grid.
     @param xcFunctional the type of exchange-correlation functional.
     @param comm the MPI communicator.
     */
    void generate(const CMolecule& molecule, const CMolecularBasis& basis, const CMolecularGrid& molGrid, const xcfun xcFunctional, MPI_Comm comm);
};

#endif /* DensityGridDriver_hpp */
