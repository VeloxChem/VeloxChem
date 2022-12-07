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

#include <mpi.h>

#include <cstdint>
#include <vector>

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
     Distributes spin-restriced density values into density grid.

     @param densityGrid the pointer to density grid object.
     @param aoDensityMatrix the AO density matrix.
     @param gtoValues the GTOs values buffer.
     @param gridOffset the offset of grid points batch in molecular grid.
     @param gridBlockPosition the position of grid block in GTOs values grid.
     @param nGridPoints the number of grid points in grid points batch.
     */
    void _distRestrictedDensityValuesForLda(CDensityGrid*              densityGrid,
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
    void _distRestrictedDensityValuesForGga(CDensityGrid*              densityGrid,
                                            const CAODensityMatrix*    aoDensityMatrix,
                                            const CMemBlock2D<double>& gtoValues,
                                            const CMemBlock2D<double>& gtoValuesX,
                                            const CMemBlock2D<double>& gtoValuesY,
                                            const CMemBlock2D<double>& gtoValuesZ,
                                            const int32_t              gridOffset,
                                            const int32_t              gridBlockPosition,
                                            const int32_t              nGridPoints) const;

    /**
     Creates density grid for pair-density LDA case.

     @param densityGrid the density grid object.
     @param aoDensityMatrix the AO density matrices.
     @param twoDM the "active" MO two-body density matrix.
     @param activeMOs the active MO coefficients.
     @param nActive the number of active orbitals.
     @param gtoContainer TODO
     @param gridCoordinatesX TODO
     @param gridCoordinatesY TODO
     @param gridCoordinatesZ TODO
     @param gridOffset the offset of grid points batch in molecular grid.
     @param nGridPoints the number of grid points in grid points batch.
     */
    void _genBatchOfPairDensityGridPointsForLda(CDensityGrid*           densityGrid,
                                                const CAODensityMatrix* aoDensityMatrix,
                                                const double*           twoDM,
                                                const double*           activeMOs,
                                                const int32_t           nActive,
                                                const CGtoContainer*    gtoContainer,
                                                const double*           gridCoordinatesX,
                                                const double*           gridCoordinatesY,
                                                const double*           gridCoordinatesZ,
                                                const int32_t           gridOffset,
                                                const int32_t           nGridPoints) const;

    /**
     Creates density grid for pair-density GGA case.

     @param densityGrid the density grid object.
     @param aoDensityMatrix the AO density matrices.
     @param twoDM the "active" MO two-body density matrix.
     @param activeMOs the active MO coefficients.
     @param nActive the number of active orbitals.
     @param gtoContainer TODO
     @param gridCoordinatesX TODO
     @param gridCoordinatesY TODO
     @param gridCoordinatesZ TODO
     @param gridOffset the offset of grid points batch in molecular grid.
     @param nGridPoints the number of grid points in grid points batch.
     */
    void _genBatchOfPairDensityGridPointsForGga(CDensityGrid*           densityGrid,
                                                const CAODensityMatrix* aoDensityMatrix,
                                                const double*           twoDM,
                                                const double*           activeMOs,
                                                const int32_t           nActive,
                                                const CGtoContainer*    gtoContainer,
                                                const double*           gridCoordinatesX,
                                                const double*           gridCoordinatesY,
                                                const double*           gridCoordinatesZ,
                                                const int32_t           gridOffset,
                                                const int32_t           nGridPoints) const;

    /**
     Distributes density and pair-density values into density grid
     or translate them into artificial alpha and beta densities.

     @param densityGrid the pointer to density grid object.
     @param twoDM the "active" MO two-body density matrix.
     @param activeMOs the active MO coefficients.
     @param nActive the number of active orbitals.
     @param nAOs TODO
     @param gtoValues the GTOs values buffer.
     @param gridOffset the offset of grid points batch in molecular grid.
     @param gridBlockPosition the position of grid block in GTOs values grid.
     @param nGridPoints the number of grid points in grid points batch.
     */
    void _distPairDensityValuesForLda(CDensityGrid*              densityGrid,
                                      const double*              twoDM,
                                      const double*              activeMOs,
                                      const int32_t              nActive,
                                      const int32_t              nAOs,
                                      const CMemBlock2D<double>& gtoValues,
                                      const int32_t              gridOffset,
                                      const int32_t              gridBlockPosition,
                                      const int32_t              nGridPoints) const;

    /**
     Distributes density and pair-density values into density grid
     or translate them into artificial alpha and beta densities.

     @param densityGrid the pointer to density grid object.
     @param twoDM the "active" MO two-body density matrix.
     @param activeMOs the active MO coefficients.
     @param nActive the number of active orbitals.
     @param nAOs TODO
     @param gtoValues the GTOs values buffer.
     @param gtoValuesX the GTOs gradient along X axis values buffer.
     @param gtoValuesY the GTOs gradient along Y axis values buffer.
     @param gtoValuesZ the GTOs gradient along Z axis values buffer.
     @param gridOffset the offset of grid points batch in molecular grid.
     @param gridBlockPosition the position of grid block in GTOs values grid.
     @param nGridPoints the number of grid points in grid points batch.
     */
    void _distPairDensityValuesForGga(CDensityGrid*              densityGrid,
                                      const double*              twoDM,
                                      const double*              activeMOs,
                                      const int32_t              nActive,
                                      const int32_t              nAOs,
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
     Generates partitioned density and on-top-pair density grid for given molecule
     and type of exchange-correlation functional. Density grid generation is distributed
     within domain of MPI communicator.

     @param aoDensityMatrix the AO density matrix.
     @param twoDM the "active" MO two-body density matrix.
     @param activeMOs the active MO coefficients.
     @param nActive the number of active orbitals.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param molecularGrid the molecular grid.
     @param xcFunctional the type of exchange-correlation functional.
     @return the density grid object.
     */
    CDensityGrid generatePdftGrid(const CAODensityMatrix& aoDensityMatrix,
                                  const double*           twoDM,
                                  const double*           activeMOs,
                                  const int32_t           nActive,
                                  const CMolecule&        molecule,
                                  const CMolecularBasis&  basis,
                                  const CMolecularGrid&   molecularGrid,
                                  const xcfun             xcFunctional);
};

#endif /* DensityGridDriver_hpp */
