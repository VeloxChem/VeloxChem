//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef GridDriver_hpp
#define GridDriver_hpp

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
     */
    CGridDriver();

    /**
     Sets accuracy level for grid generation. Level: 1-8, where 1 is coarse
     grid, 5 is ultrafine grid, 8 special benchmarking grid.

     @param gridLevel the accuracy level of generated grid.
     */
    auto setLevel(const int gridLevel) -> void;

    /**
     Generates molecular grid for molecule.

     @param molecule the molecule.
     @param rank the MPI rank.
     @param nnodes the number of MPI processes.
     @return the molecular grid object.
     */
    auto generate_local_grid(const CMolecule& molecule, const int rank, const int nnodes) const -> CMolecularGrid;
};

#endif /* GridDriver_hpp */
