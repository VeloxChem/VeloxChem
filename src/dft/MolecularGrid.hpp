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

#ifndef MolecularGrid_hpp
#define MolecularGrid_hpp

#include <cstdint>
#include <ostream>
#include <array>

#include <mpi.h>

#include "MemBlock2D.hpp"

/**
 Class CMolecularGrid class generates molecular grid.

 @author Z. Rinkevicius
 */
class CMolecularGrid
{
    /**
     The distribution status of molecular grid object.
     */
    bool _isDistributed;

    /**
     The grid points (coordinates, weights).
     */
    CMemBlock2D<double> _gridPoints;

   public:
    /**
     Creates an empty molecular grid object.
     */
    CMolecularGrid();

    /**
     Creates a molecular grid object.

     @param gridPoints the 2D memory block object with grid points data.
     */
    CMolecularGrid(const CMemBlock2D<double>& gridPoints);

    /**
     Creates a molecular grid object by copying other molecular grid object.

     @param source the molecular grid object.
     */
    CMolecularGrid(const CMolecularGrid& source);

    /**
     Creates a molecular grid object by by moving other molecular grid object.

     @param source the molecular grid object.
     */
    CMolecularGrid(CMolecularGrid&& source) noexcept;

    /**
     Destroys a molecular grid object.
     */
    ~CMolecularGrid();

    /**
     Assigns a molecular grid object by copying other molecular grid object.

     @param source the molecular grid object.
     */
    CMolecularGrid& operator=(const CMolecularGrid& source);

    /**
     Assigns a molecular grid object by moving other molecular grid object.

     @param source the molecular grid object.
     */
    CMolecularGrid& operator=(CMolecularGrid&& source) noexcept;

    /**
     Compares molecular grid object with other molecular grid object.

     @param other the molecular grid object.
     @return true if molecular grid objects are equal, false otherwise.
     */
    bool operator==(const CMolecularGrid& other) const;

    /**
     Compares molecular grid object with other molecular grid object.

     @param other the molecular grid object.
     @return true if molecular grid objects are not equal, false otherwise.
     */
    bool operator!=(const CMolecularGrid& other) const;
    
    
    /**
     Reduces size of molecular grid by slicing all grid points beoynd given number of grid points.

     @param nGridPoints the number of grid points.
     */
    void slice(const int32_t nGridPoints);

    /**
     Gets number of grid points in molecular grid object.

     @return the number of grid points.
     */
    int32_t getNumberOfGridPoints() const;

    /**
     Gets Cartesian X coordinates of grid points in molecular grid object.

     @return the constant pointer to Cartesian X coordinates of grid points.
     */
    const double* getCoordinatesX() const;
    
    /**
     Gets Cartesian X coordinates of grid points in molecular grid object.
     
     @return the pointer to Cartesian X coordinates of grid points.
     */
    double* getCoordinatesX();

    /**
     Gets Cartesian Y coordinates of grid points in molecular grid object.

     @return the constant  pointer to Cartesian Y coordinates of grid points.
     */
    const double* getCoordinatesY() const;
    
    /**
     Gets Cartesian Y coordinates of grid points in molecular grid object.
     
     @return the pointer to Cartesian Y coordinates of grid points.
     */
    double* getCoordinatesY();

    /**
     Gets Cartesian Z coordinates of grid points in molecular grid object.

     @return the constant pointer to Cartesian Z coordinates of grid points.
     */
    const double* getCoordinatesZ() const;
    
    /**
     Gets Cartesian Z coordinates of grid points in molecular grid object.
     
     @return the pointer to Cartesian Z coordinates of grid points.
     */
    double* getCoordinatesZ();

    /**
     Gets weights of grid points in molecular grid object.

     @return the constant pointer to weights of grid points.
     */
    const double* getWeights() const;
    
    /**
     Gets weights of grid points in molecular grid object.
     
     @return the pointer to weights of grid points.
     */
    double* getWeights();

    /**
     Distributes grid points data across a molecular grid objects associated
     with MPI processes within domain of MPI communacator and sets distribution
     flag to true.

     @param rank the rank of MPI process.
     @param nodes the number of nodes in MPI domain.
     @param comm he MPI communicator.
     */
    void distribute(int32_t rank, int32_t nodes, MPI_Comm comm);

    /**
     Broadcasts grid points data across a molecular grid objects within domain
     of MPI communacator.

     @param rank the rank of MPI process.
     @param comm he MPI communicator.
     */
    void broadcast(int32_t rank, MPI_Comm comm);
    
    /**
     Reads raw blocked grid data from DALTON program quad file. NOTE: maximum allowed number of grid points
     is 20M.

     @param fileName the name of binary quad file.
     */
    void read_blocked_grid(const std::string& fileName);
    
    /**
     Reads raw blocked grid data from text file in format: (grid point index, rx, ry, rz, w). NOTE: maximum allowed
     number of grid points is 20M.
     
     @param fileName the name of text file.
     */
    void read_raw_grid(const std::string& fileName);
    
    /**
     Writes raw grid data to text file in format: (grid point index, rx, ry, rz, w).

     @param fileName the name of text file.
     */
    void write_raw_grid(const std::string fileName) const;
    
    /**
     Gets spatial extent of molecular grid.

     @return the spatial extent (min x, min y, min z, max x, max y, max z).
     */
    std::array<double, 6> getSpatialExtent() const;

    /**
     Converts molecular grid object to text and insert it into output text
     stream.

     @param output the output text stream.
     @param source the molecular grid.
     */
    friend std::ostream& operator<<(std::ostream& output, const CMolecularGrid& source);
};

#endif /* MolecularGrid_hpp */
