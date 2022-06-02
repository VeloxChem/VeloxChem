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

#ifndef DensityGradientGridDriver_hpp
#define DensityGradientGridDriver_hpp

#include <mpi.h>

#include <cstdint>
#include <vector>

#include "AODensityMatrix.hpp"
#include "DensityGrid.hpp"
#include "GtoContainer.hpp"
#include "MolecularBasis.hpp"
#include "MolecularGrid.hpp"
#include "Molecule.hpp"

/**
 Class CDensityGradientGridDriver generates density nuclear gradient grid for usage in numerical
 integration.

 @author Z. Rinkevicius
 */
class CDensityGradientGridDriver
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
     Creates density grid for spin-restricted LDA case.

     @param densityGrid the density grid object.
     @param aoDensityMatrix the AO density matrices.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param molecularGrid the distributed molecular grid.
     @param iAtom the index of selected atom.
     */
    void _genRestrictedDensityForLda(CDensityGrid&           densityGrid,
                                     const CAODensityMatrix& aoDensityMatrix,
                                     const CMolecule&        molecule,
                                     const CMolecularBasis&  basis,
                                     const CMolecularGrid&   molecularGrid,
                                     const int32_t           iAtom) const;

    /**
     Generates batch of spin restricted density grid points for LDA case.

     @param densityGrid the pointer to density grid object.
     @param aoDensityMatrix the AO density matrix.
     @param gtoContainer the GTOs container.
     @param atmGtoContainer the GTOs container for specific atom
     @param gridCoordinatesX the vector of Cartesian X coordinates of grid
     points.
     @param gridCoordinatesY the vector of Cartesian Y coordinates of grid
     points.
     @param gridCoordinatesZ the vector of Cartesian Y coordinates of grid
     points.
     @param gridOffset the batch offset in vector grid points.
     @param nGridPoints the number of grid points in batch.
     */
    void _genBatchOfRestrictedDensityGridPointsForLda(CDensityGrid*           densityGrid,
                                                      const CAODensityMatrix* aoDensityMatrix,
                                                      const CGtoContainer*    gtoContainer,
                                                      const CGtoContainer*    atmGtoContainer,
                                                      const double*           gridCoordinatesX,
                                                      const double*           gridCoordinatesY,
                                                      const double*           gridCoordinatesZ,
                                                      const int32_t           gridOffset,
                                                      const int32_t           nGridPoints) const;

    /**
     Distributes spin-restriced density values into density grid.

     @param densityGrid the pointer to density grid object.
     @param aoDensityMatrix the AO density matrix.
     @param braGtoValues the GTOs values buffer on bra side.
     @param ketGtoValuesX the GTOs gradient along X axis values buffer on ket side.
     @param ketGtoValuesY the GTOs gradient along Y axis values buffer on ket side.
     @param ketGtoValuesZ the GTOs gradient along Z axis values buffer on ket side.
     @param gridOffset the offset of grid points batch in molecular grid.
     @param gridBlockPosition the position of grid block in GTOs values grid.
     @param nGridPoints the number of grid points in grid points batch.
     */
    void _distRestrictedDensityValuesForLda(CDensityGrid*              densityGrid,
                                            const CAODensityMatrix*    aoDensityMatrix,
                                            const CMemBlock2D<double>& braGtoValues,
                                            const CMemBlock2D<double>& ketGtoValuesX,
                                            const CMemBlock2D<double>& ketGtoValuesY,
                                            const CMemBlock2D<double>& ketGtoValuesZ,
                                            const int32_t              gridOffset,
                                            const int32_t              gridBlockPosition,
                                            const int32_t              nGridPoints) const;

    /**
     Creates density grid for spin-restricted GGA case.

     @param densityGrid the density grid object.
     @param aoDensityMatrix the AO density matrices.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param molecularGrid the distributed molecular grid.
     @param iAtom the index of selected atom.
     */
    void _genRestrictedDensityForGga(CDensityGrid&           densityGrid,
                                     const CAODensityMatrix& aoDensityMatrix,
                                     const CMolecule&        molecule,
                                     const CMolecularBasis&  basis,
                                     const CMolecularGrid&   molecularGrid,
                                     const int32_t           iAtom) const;

    /**
     Generates batch of spin restricted density grid points for GGA case.

     @param densityGrid the pointer to density grid object.
     @param aoDensityMatrix the AO density matrix.
     @param gtoContainer the GTOs container.
     @param atmGtoContainer the GTOs container for specific atom
     @param gridCoordinatesX the vector of Cartesian X coordinates of grid
     points.
     @param gridCoordinatesY the vector of Cartesian Y coordinates of grid
     points.
     @param gridCoordinatesZ the vector of Cartesian Y coordinates of grid
     points.
     @param gridOffset the batch offset in vector grid points.
     @param nGridPoints the number of grid points in batch.
     */

    void _genBatchOfRestrictedDensityGridPointsForGga(CDensityGrid*           densityGrid,
                                                      const CAODensityMatrix* aoDensityMatrix,
                                                      const CGtoContainer*    gtoContainer,
                                                      const CGtoContainer*    atmGtoContainer,
                                                      const double*           gridCoordinatesX,
                                                      const double*           gridCoordinatesY,
                                                      const double*           gridCoordinatesZ,
                                                      const int32_t           gridOffset,
                                                      const int32_t           nGridPoints) const;

    /**
     Distributes spin-restriced density values into density grid.

     @param densityGrid the pointer to density grid object.
     @param aoDensityMatrix the AO density matrix.
     @param braGtoValues the GTOs values buffer on bra side.
     @param braGtoValuesX the GTOs gradient along X axis values buffer on bra side.
     @param braGtoValuesY the GTOs gradient along Y axis values buffer on bra side.
     @param braGtoValuesZ the GTOs gradient along Z axis values buffer on bra side.
     @param ketGtoValuesX the GTOs gradient along X axis values buffer on ket side.
     @param ketGtoValuesY the GTOs gradient along Y axis values buffer on ket side.
     @param ketGtoValuesZ the GTOs gradient along Z axis values buffer on ket side.
     @param gridOffset the offset of grid points batch in molecular grid.
     @param gridBlockPosition the position of grid block in GTOs values grid.
     @param nGridPoints the number of grid points in grid points batch.
     */

    void _distRestrictedDensityValuesForGga(CDensityGrid*              densityGrid,
                                            const CAODensityMatrix*    aoDensityMatrix,
                                            const CMemBlock2D<double>& braGtoValues,
                                            const CMemBlock2D<double>& braGtoValuesX,
                                            const CMemBlock2D<double>& braGtoValuesY,
                                            const CMemBlock2D<double>& braGtoValuesZ,
                                            const CMemBlock2D<double>& ketGtoValuesX,
                                            const CMemBlock2D<double>& ketGtoValuesY,
                                            const CMemBlock2D<double>& ketGtoValuesZ,
                                            const CMemBlock2D<double>& ketGtoValuesXX,
                                            const CMemBlock2D<double>& ketGtoValuesXY,
                                            const CMemBlock2D<double>& ketGtoValuesXZ,
                                            const CMemBlock2D<double>& ketGtoValuesYY,
                                            const CMemBlock2D<double>& ketGtoValuesYZ,
                                            const CMemBlock2D<double>& ketGtoValuesZZ,
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
     Creates a density gradient grid driver object using MPI info.

     @param comm the MPI communicator.
     */
    CDensityGradientGridDriver(MPI_Comm comm);

    /**
     Destroys a density gradient grid driver object.
     */
    ~CDensityGradientGridDriver();

    /**
     Generates partitioned density gradient grid for given atom in molecule.

     @param aoDensityMatrix the AO density matrix.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param molecularGrid the molecular grid.
     @param xcFunctional the exchange-correlation functional type.
     @param iAtom the index of atom in molecule.
     @return the density gradeint grid object.
     */
    CDensityGrid generate(const CAODensityMatrix& aoDensityMatrix,
                          const CMolecule&        molecule,
                          const CMolecularBasis&  basis,
                          const CMolecularGrid&   molecularGrid,
                          const xcfun             xcFunctional,
                          const int32_t           iAtom);
};

#endif /* DensityGradientGridDriver_hpp */
