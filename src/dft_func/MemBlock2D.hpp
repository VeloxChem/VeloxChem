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

#ifndef MemBlock2D_hpp
#define MemBlock2D_hpp

#include <vector>

/**
 Templated class CMemBlock2D manages 2D memory block allocation, manipulation,
 and deallocation.
 */
class CMemBlock2D
{
    /**
     The contiguous memory block.
     */
    std::vector<double> _data;

    /**
     The original sizes of data chunks in memory block.
     */
    std::vector<int> _originalSizes;

    /**
     The padded sizes of data chunks in memory block.
     */
    std::vector<int> _paddedSizes;

    /**
     The start positions of data chunks in memory block.
     */
    std::vector<int> _positions;

    /**
     The number of elements in memory block.
     */
    int _nElements;

    /**
     Sets uniform original size for all data chunks.

     @param nElements the number of elements in data chunk.
     @param nBlocks the number of data chunks.
     */
    void _setOriginalSizes(const int nElements, const int nBlocks);

    /**
     Sets dimensions, padding for all data chunks.
     */
    void _setDimensions();

   public:
    /**
     Creates an empty 2D memory block object.
     */
    CMemBlock2D();

    /**
     Creates an 2D memory block object.

     @param nElements the number of elements in data chunk.
     @param nBlocks the number of data chunks.
     */
    CMemBlock2D(const int nElements, const int nBlocks);

    /**
     Creates a 2D memory block object by copying other 2D memory block object.

     @param source the 2D memory block object.
     */
    CMemBlock2D(const CMemBlock2D& source);

    /**
     Creates a 2D memory block object by moving other 2D memory block object.

     @param source the 2D memory block object.
     */
    CMemBlock2D(CMemBlock2D&& source) noexcept;

    /**
     Destroys a 2D memory block object.
     */
    ~CMemBlock2D();

    /**
     Assigns a 2D memory block object by copying other 2D memory block object.

     @param source the 2D memory block object.
     */
    CMemBlock2D& operator=(const CMemBlock2D& source);

    /**
     Assigns a 2D memory block object by moving other 2D memory block object.

     @param source the 2D memory block object.
     */
    CMemBlock2D& operator=(CMemBlock2D&& source) noexcept;

    /**
     Compares 2D memory block object with other 2D memory block object.

     @param other the 2D memory block object.
     @return true if 2D memory block objects are equal, false otherwise.
     */
    bool operator==(const CMemBlock2D& other) const;

    /**
     Compares 2D memory block object with other 2D memory block object.

     @param other the 2D memory block object.
     @return true if 2D memory block objects are not equal, false otherwise.
     */
    bool operator!=(const CMemBlock2D& other) const;

    /**
     Sets all elements of contiguous memory block to zero.
     */
    void zero();

    /**
     Checks if 2D memory block object is empty.

     @return true if 2D memory block object is empty, false otherwise.
     */
    bool isEmpty() const;

    /**
     Gets pointer to specific data chunk in 2D memory block.

     @param iBlock the index of data chunk.
     @return the pointer to data chunk.
     */
    double* data(const int iBlock);

    /**
     Gets constant pointer to specific data chunk in 2D memory block.

     @param iBlock the index of data chunk.
     @return the pointer to data chunk.
     */
    const double* data(const int iBlock) const;

    /**
     Gets pointer to specific element of selected data chunk in 2D memory block.

     @param iBlock the index of data chunk.
     @param iElement the index of element in memory block.
     @return the pointer to data chunk.
     */
    double* data(const int iBlock, const int iElement);

    /**
     Gets constant pointer to specific elment of selected data chunk in 2D
     memory block.

     @param iBlock the index of data chunk.
     @param iElement the index of element in memory block.
     @return the pointer to data chunk.
     */
    const double* data(const int iBlock, const int iElement) const;

    /**
     Gets number of elements in specific data chunk.

     @param iBlock the index of data chunk.
     @return the number of elements in data chunk.
     */
    int size(const int iBlock) const;

    /**
     Gets number of elements in specific pitched data chunk.

     @param iBlock the index of data chunk.
     @return the number of elements in pitched data chunk.
     */
    int pitched_size(const int iBlock) const;

    /**
     Gets number of data chunk.

     @return the number of data chunk.
     */
    int blocks() const;
};

#endif /* MemBlock2D_hpp */
