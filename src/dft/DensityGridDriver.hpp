//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef DensityGridDriver_hpp
#define DensityGridDriver_hpp

#include <cstdint>
#include <vector>

#include "mpi.h"

#include "Molecule.hpp"
#include "MolecularBasis.hpp"
#include "MolecularGrid.hpp"
#include "XCFuncType.hpp"
#include "ExecMode.hpp"
#include "GtoContainer.hpp"
#include "VecMemBlocks.hpp"

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
     @param oStream the output stream.
     @param comm the MPI communicator.
     */
    void _genDensityGridOnCPU(const CMolecule&       molecule,
                              const CMolecularBasis& basis,
                              const CMolecularGrid&  molGrid,
                              const xcfun            xcFunctional,
                                    COutputStream&   oStream,
                                    MPI_Comm         comm);
    
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
     */
    void _genBatchOfDensityGridPoints(const CGtoContainer* gtoContainer,
                                      const double*        gridCoordinatesX,
                                      const double*        gridCoordinatesY,
                                      const double*        gridCoordinatesZ,
                                      const int32_t        gridOffset,
                                      const int32_t        nGridPoints);
    
    /**
     Computes screening factors  vector for primitive Gaussian functions for
     specific grid point using GTOs container.

     @param gtoContainer the GTOs container.
     @param gridCoordinateX the Cartesian X coordinate of grid point.
     @param gridCoordinateY the Cartesian Y coordinate of grid point.
     @param gridCoordinateZ the Cartesian Z coordinate of grid point.
     @param screenFactors the screening factors vector for each GTOs block
            object in GTOs container.
     */
    void _compScreeningFactors(const CGtoContainer&        gtoContainer,
                               const double                gridCoordinateX,
                               const double                gridCoordinateY,
                               const double                gridCoordinateZ,
                                     CVecMemBlock<double>& screenFactors);
    
    /**
     Gets scaling factor of screening factor for specific angular momemntum.

     @param angularMomentum the angular momentum.
     @param maxRadius the maximum radius.
     @return the scaling factor.
     */
    double _getScaleFactor(const int32_t angularMomentum,
                           const double  maxRadius) const;
    
public:
    
    /**
     Creates a density grid driver object using MPI info.
     
     @param globRank the the rank of MPI process.
     @param globNodes the total number of MPI processes.
     @param runMode the execution mode.
     @param comm the MPI communicator.
     */
    CDensityGridDriver(const int32_t  globRank,
                       const int32_t  globNodes,
                       const execmode runMode,
                             MPI_Comm comm);
    
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
     @param oStream the output stream.
     @param comm the MPI communicator.
     */
    void generate(const CMolecule&       molecule,
                  const CMolecularBasis& basis,
                  const CMolecularGrid&  molGrid,
                  const xcfun            xcFunctional,
                        COutputStream&   oStream,
                        MPI_Comm         comm);
};

#endif /* DensityGridDriver_hpp */
