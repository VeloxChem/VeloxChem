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

#ifndef DensityGridData2D_hpp
#define DensityGridData2D_hpp

#include <vector>

/**
 Class CDensityGridData2D manages 2D memory block for density grid.
 */
class CDensityGridData2D
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
    CDensityGridData2D();

    /**
     Creates an 2D memory block object.

     @param nElements the number of elements in data chunk.
     @param nBlocks the number of data chunks.
     */
    CDensityGridData2D(const int nElements, const int nBlocks);

    /**
     Creates a 2D memory block object by copying other 2D memory block object.

     @param source the 2D memory block object.
     */
    CDensityGridData2D(const CDensityGridData2D& source);

    /**
     Creates a 2D memory block object by moving other 2D memory block object.

     @param source the 2D memory block object.
     */
    CDensityGridData2D(CDensityGridData2D&& source) noexcept;

    /**
     Destroys a 2D memory block object.
     */
    ~CDensityGridData2D();

    /**
     Assigns a 2D memory block object by copying other 2D memory block object.

     @param source the 2D memory block object.
     */
    CDensityGridData2D& operator=(const CDensityGridData2D& source);

    /**
     Assigns a 2D memory block object by moving other 2D memory block object.

     @param source the 2D memory block object.
     */
    CDensityGridData2D& operator=(CDensityGridData2D&& source) noexcept;

    /**
     Compares 2D memory block object with other 2D memory block object.

     @param other the 2D memory block object.
     @return true if 2D memory block objects are equal, false otherwise.
     */
    bool operator==(const CDensityGridData2D& other) const;

    /**
     Compares 2D memory block object with other 2D memory block object.

     @param other the 2D memory block object.
     @return true if 2D memory block objects are not equal, false otherwise.
     */
    bool operator!=(const CDensityGridData2D& other) const;

    /**
     Sets all elements of contiguous memory block to zero.
     */
    void zero();

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

#endif /* DensityGridData2D_hpp */
