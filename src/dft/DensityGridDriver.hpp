//
//                           VELOXCHEM 1.0-RC2
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

#ifndef DensityGridDriver_hpp
#define DensityGridDriver_hpp

#include <cstdint>
#include <vector>

#include <mpi.h>

#include "ExecMode.hpp"
#include "XCFuncType.hpp"

class CAODensityMatrix;
class CDensityGrid;
class CGtoContainer;
template <typename T>
class CMemBlock2D;
class CMolecularBasis;
class CMolecularGrid;
class CMolecule;
class CSphericalMomentum;

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
    void _genRestrictedDensityForLda(      CDensityGrid&     densityGrid,
                                     const CAODensityMatrix& aoDensityMatrix,
                                     const CMolecule&        molecule,
                                     const CMolecularBasis&  basis,
                                     const CMolecularGrid&   molecularGrid) const;
    

    /**
     Creates density grid for spin-unrestricted LDA case.
     
     @param densityGrid the density grid object.
     @param aoDensityMatrix the AO density matrices.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param molecularGrid the distributed molecular grid.
     */
    void _genUnrestrictedDensityForLda(      CDensityGrid&     densityGrid,
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
    void _genRestrictedDensityForGga(      CDensityGrid&     densityGrid,
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
    void _genUnrestrictedDensityForGga(      CDensityGrid&     densityGrid,
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
    void _genBatchOfRestrictedDensityGridPointsForLda(      CDensityGrid*     densityGrid,
                                                      const CAODensityMatrix* aoDensityMatrix,
                                                      const CGtoContainer*    gtoContainer,
                                                      const double*           gridCoordinatesX,
                                                      const double*           gridCoordinatesY,
                                                      const double*           gridCoordinatesZ,
                                                      const int32_t           gridOffset,
                                                      const int32_t           nGridPoints) const;
    

  /**
      Generates batch of spin-unrestricted density grid points for LDA case.

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
    void _genBatchOfUnrestrictedDensityGridPointsForLda(       CDensityGrid*     densityGrid,
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
    void _genBatchOfRestrictedDensityGridPointsForGga(      CDensityGrid*     densityGrid,
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
    void _genBatchOfUnrestrictedDensityGridPointsForGga(      CDensityGrid*     densityGrid,
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
    void _distRestrictedDensityValuesForLda(      CDensityGrid*        densityGrid,
                                            const CAODensityMatrix*    aoDensityMatrix,
                                            const CMemBlock2D<double>& gtoValues,
                                            const int32_t              gridOffset,
                                            const int32_t              gridBlockPosition,
                                            const int32_t              nGridPoints) const;

    /**
     Distributes spin-unrestriced density values into density grid.

     @param densityGrid the pointer to density grid object.
     @param aoDensityMatrix the AO density matrix.
     @param gtoValues the GTOs values buffer.
     @param gridOffset the offset of grid points batch in molecular grid.
     @param gridBlockPosition the position of grid block in GTOs values grid.
     @param nGridPoints the number of grid points in grid points batch.
     */
    void _distUnrestrictedDensityValuesForLda(      CDensityGrid*        densityGrid,
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
    void _distRestrictedDensityValuesForGga(      CDensityGrid*        densityGrid,
                                            const CAODensityMatrix*    aoDensityMatrix,
                                            const CMemBlock2D<double>& gtoValues,
                                            const CMemBlock2D<double>& gtoValuesX,
                                            const CMemBlock2D<double>& gtoValuesY,
                                            const CMemBlock2D<double>& gtoValuesZ,
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
    void _distUnrestrictedDensityValuesForGga(      CDensityGrid*        densityGrid,
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


    CDensityGrid pdft(const CAODensityMatrix& aoDensityMatrix,
                         double* twoDM,
                         double* activeMOs,
                         int nActive,
                             const CMolecule&        molecule,
                             const CMolecularBasis&  basis,
                             const CMolecularGrid&   molecularGrid,
                             const xcfun             xcFunctional);

    void _PDFT_Lda(      CDensityGrid*     densityGrid,
                         const CAODensityMatrix* aoDensityMatrix,
                         double* twoDM,
                         double* activeMOs,
                         int nActive,
                         const CGtoContainer*    gtoContainer,
                         const double*           gridCoordinatesX,
                         const double*           gridCoordinatesY,
                         const double*           gridCoordinatesZ,
                         const int32_t           gridOffset,
                         const int32_t           nGridPoints) const;

    void _distPDFT_LDA(      CDensityGrid*        densityGrid,
                                        double* twoDM,
                                        double* activeMOs,
                                        int nActive,
                                        int nAOs,
                                        const CMemBlock2D<double>& gtoValues,
                                        const int32_t              gridOffset,
                                        const int32_t              gridBlockPosition,
                                        const int32_t              nGridPoints) const;

};

#endif /* DensityGridDriver_hpp */
