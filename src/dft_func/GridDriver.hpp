//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
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

#ifndef GridDriver_hpp
#define GridDriver_hpp

#include <mpi.h>

#include <cstdint>

#include "MolecularGrid.hpp"
#include "Molecule.hpp"
#include "Point.hpp"

/**
 Class CGridDriver generates grid points data for usage in numerical
 integration.
 */
class CGridDriver
{
    /**
     The accuracy level of grid (from 1 to 6).
     */
    int _gridLevel;

    /**
     The threshold of weights screening.
     */
    double _thresholdOfWeight;

    /**
     The MPI communicator.
     */
    MPI_Comm _locComm;

    /**
     Determines number of radial grid points for specific chemical element.

     @param idElemental the chemical element number.
     @return the number of radial points.
     */
    auto _getNumberOfRadialPoints(const int idElemental) const -> int;

    /**
     Determines number of angular grid points for specific chemical element.

     @param idElemental the chemical element number.
     @return the number of angular points.
     */
    auto _getNumberOfAngularPoints(const int idElemental) const -> int;

    /**
     Creates molecular grid on master node by generating fraction of grid
     points on each MPI process within domain of MPI communicator. Grid points
     are generated using only CPUs.

     @param molecule the molecule.
     @return the molecular grid object.
     */
    auto _genGridPoints(const CMolecule& molecule) const -> CMolecularGrid;

    /**
     Gets size of grid points batch.

     @param idsElemental the vector of chemical elements identifiers.
     @param offset the  in vector of chemical elements identifiers.
     @param nAtoms the number of atoms in batch.
     @return the number of grid points.
     */
    auto _getBatchSize(const int* idsElemental, const int offset, const int nAtoms) const -> int;

    /**
     Generates grid points for specific atom in molecule.

     @param rawGridPoints the raw grid points.
     @param minDistance the distance to closest neighbouring atom.
     @param gridOffset the atom grid points offset in raw grid points.
     @param atomCoordinates the Cartesian coordinates of atoms.
     @param nAtoms the number of atoms.
     @param idElemental the chemical element identifier of atom.
     @param idAtomic the index of atom.
     */
    auto _genAtomGridPoints(CDenseMatrix*   rawGridPoints,
                            const double    minDistance,
                            const int   gridOffset,
                            const TPoint<double>* atomCoordinates,
                            const int   nAtoms,
                            const int   idElemental,
                            const int   idAtomic) const -> void;

    /**
     Prunes raw grid points by screening weights and discarding all grid point
     with weights bellow cutoff threshold.

     @param rawGridPoints the raw grid points.
     @return the number of pruned grid points.
     */
    auto _screenRawGridPoints(CDenseMatrix* rawGridPoints) const -> int;

   public:
    /**
     Creates a grid driver object.

     @param comm the MPI communicator.
     */
    CGridDriver(MPI_Comm comm);

    /**
     Sets accuracy level for grid generation. Level: 1-8, where 1 is coarse
     grid, 5 is ultrafine grid, 8 special benchmarking grid.

     @param gridLevel the accuracy level of generated grid.
     */
    auto setLevel(const int gridLevel) -> void;

    /**
     Generates molecular grid for molecule. Errors are printed to output stream.
     Grid generation is distributed within domain of MPI communicator.

     @param molecule the molecule.
     @return the molecular grid object.
     */
    auto generate(const CMolecule& molecule) const -> CMolecularGrid;
};

#endif /* GridDriver_hpp */
