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

     @param denGrid the density grid object. 
     @param density the AO density matrices.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param molGrid the distributed molecular grid.
     @param xcFunctional the exchange-correlation functional type.
     */
    void _genDensityGridOnCPU(      CDensityGrid&     denGrid,
                              const CAODensityMatrix& density,
                              const CMolecule&        molecule,
                              const CMolecularBasis&  basis,
                              const CMolecularGrid&   molGrid,
                              const xcfun             xcFunctional);

    /**
     Generates batch of density grid points.

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
     @param xcFunctional the exchange-correlation functional type.
     */
    void _genBatchOfDensityGridPoints(      CDensityGrid*     densityGrid,
                                      const CAODensityMatrix* aoDensityMatrix,
                                      const CGtoContainer*    gtoContainer,
                                      const double*           gridCoordinatesX,
                                      const double*           gridCoordinatesY,
                                      const double*           gridCoordinatesZ,
                                      const int32_t           gridOffset,
                                      const int32_t           nGridPoints,
                                      const xcfun             xcFunctional);
    
    
    /**
     Computes density at grid point for pair of GTO blocks.
     
     @param densityGrid the pointer to density grid object.
     @param aoDensityMatrix the AO density matrix.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
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
    void _compDensityForGtoBlocks(      CDensityGrid*     densityGrid,
                                  const CAODensityMatrix* aoDensityMatrix,
                                  const CGtoBlock&        braGtoBlock,
                                  const CGtoBlock&        ketGtoBlock,
                                  const double*           gridCoordinatesX,
                                  const double*           gridCoordinatesY,
                                  const double*           gridCoordinatesZ,
                                  const int32_t           gridOffset,
                                  const int32_t           nGridPoints,
                                  const xcfun             xcFunctional);
    
    /**
     Gets number of components required for compupation of all relevant density
     contributions for specific type of exchange-correlation functions.

     @param xcFunctional the exchange-correlations functional type.
     @return the number of components.
     */
    int32_t _getNumberOfXCComponents(const xcfun xcFunctional) const;
    
    
    /**
     Computes contracted GTOs values at grid points for specific type of functional.

     @param cartGtoGridBuffer the buffer for storing primitive Cartesian GTOs values on the grid.
     @param spherGtoGridBuffer the buffer for storing contracted spherical GTOs values on the grid.
     @param gridCoordinatesX the vector of Cartesian X coordinates of grid
     points.
     @param gridCoordinatesY the vector of Cartesian Y coordinates of grid
     points.
     @param gridCoordinatesZ the vector of Cartesian Y coordinates of grid
     points.
     @param gridOffset the batch offset in vector grid points.
     @param gtoBlock the GTOs block.
     @param iContrGto the index of contracted GTO is GTOs block.
     @param xcFunctional the exchange-correlations functional type.
     */
    void _compGtoValuesOnGrid(      CMemBlock2D<double>& cartGtoGridBuffer,
                                    CMemBlock2D<double>& spherGtoGridBuffer,
                              const double*              gridCoordinatesX,
                              const double*              gridCoordinatesY,
                              const double*              gridCoordinatesZ,
                              const int32_t              gridOffset,
                              const CGtoBlock&           gtoBlock,
                              const int32_t              iContrGto,
                              const xcfun                xcFunctional) const;

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
    void _transContrGtoValues(      CMemBlock2D<double>&  spherGtoValues,
                              const CMemBlock2D<double>&  cartGtoValues,
                              const CSphericalMomentum&   spherMomentum,
                              const CMemBlock2D<int32_t>& redDimensions,
                              const int32_t               iGtoBlock,
                              const xcfun                 xcFunctional) const;
    
    
    /**
     Computes density values at given grid point for specific functional type.

     @param spherGtoValues he spherical GTOs values.
     @param gtoContainer the GTOs container.
     @param xcFunctional the exchange-correlation functional type.
     */
    void _compDensityValues(const CVecMemBlock2D<double>& spherGtoValues,
                            const CGtoContainer&          gtoContainer,
                            const xcfun                   xcFunctional) const;

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

     @param density the AO density matrix.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param molGrid the molecular grid.
     @param xcFunctional the type of exchange-correlation functional.
     @return the density grid object.
     */
    CDensityGrid generate(const CAODensityMatrix& density,
                          const CMolecule&        molecule,
                          const CMolecularBasis&  basis,
                          const CMolecularGrid&   molGrid,
                          const xcfun             xcFunctional);
};

#endif /* DensityGridDriver_hpp */
