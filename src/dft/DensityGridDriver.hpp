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

     @param spherGtoGridBuffer the buffer for storing contracted spherical GTOs values on the grid.
     @param cartGtoGridBuffer the buffer for storing primitive Cartesian GTOs values on the grid.
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
    void _compGtoValuesOnGrid(      CMemBlock2D<double>& spherGtoGridBuffer,
                                    CMemBlock2D<double>& cartGtoGridBuffer,
                              const double*              gridCoordinatesX,
                              const double*              gridCoordinatesY,
                              const double*              gridCoordinatesZ,
                              const int32_t              gridOffset,
                              const CGtoBlock&           gtoBlock,
                              const int32_t              iContrGto,
                              const xcfun                xcFunctional) const;
    
    /**
     Sets density pairs from AO density matrix.

     @param densityPairs the density pairs data.
     @param aoDensityMatrix the AO density matrix.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param isBraEqualKet the flag for equality of bra and ket GTOs blocks.
     @param iBraContrGto the index of contracted GTO on bra side.
     @param iKetContrGto the index of contracted GTO on ket side.
     @return true if density pair contains non-vanishing element, false otherwise.
     */
    bool _setDensityPair(      CMemBlock2D<double>& densityPairs,
                         const CAODensityMatrix*    aoDensityMatrix,
                         const CGtoBlock&           braGtoBlock,
                         const CGtoBlock&           ketGtoBlock,
                         const bool                 isBraEqualKet,
                         const int32_t              iBraContrGto,
                         const int32_t              iKetContrGto) const;
    
    
    /**
     Adds GTOs pair contribution to density grid for specific type of functional.

     @param densityGrid the density grid object.
     @param densityPairs the density pair data.
     @param isRestrictedDensity the flag whatever density pair data is for restricted or unrestricted AO densities.
     @param braGtoValues the GTOs values on bra side.
     @param ketGtoValues the GTOs values on ket side.
     @param gridOffset the offset og density grid.
     @param xcFunctional the exchange-correlation functional type.
     */
    void _addGtosPairContribution(      CDensityGrid*        densityGrid,
                                  const CMemBlock2D<double>& densityPairs,
                                  const bool                 isRestrictedDensity,
                                  const CMemBlock2D<double>& braGtoValues,
                                  const CMemBlock2D<double>& ketGtoValues,
                                  const int32_t              gridOffset,
                                  const xcfun                xcFunctional) const;

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
