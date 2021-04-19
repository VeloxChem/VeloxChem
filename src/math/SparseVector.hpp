//
//                           VELOXCHEM 1.0-RC
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

#ifndef SparseVector_hpp
#define SparseVector_hpp

#include <cstdint>
#include <ostream>
#include <vector>

#include "MemBlock.hpp"

/**
 Class CSparseVector stores sparse vector (zero-based indexing scheme) and
 provides set of methods for manipulating sparse vector data.

 @author Z. Rinkevicius
 */
class CSparseVector
{
    /**
     The length of vector.
     */
    int32_t _nLength;

    /**
     The values of vector elements.
     */
    CMemBlock<double> _values;

    /**
     The vector of element indexes.
     */
    CMemBlock<int32_t> _indexes;

    /**
     The maximum number of elements in values vector.
     */
    int32_t _nMaxElements;

    /**
     The current number of elements in values vector.
     */
    int32_t _nElements;

    /**
     The cut-off threshold for small elements in vector.
     */
    double _threshold;

    /**
     Sets maximum number of initial elments in sparse vector.
     */
    int32_t _setMaxNumberOfElements() const;

   public:
    /**
     Creates an empty sparse vector object.
     */
    CSparseVector();

    /**
     Creates a sparse vector object from vector of values and vector of indexes.

     @param values the vector of matrix elements.
     @param indexes the vector of indexes.
     @param nLength the length of vector.
     @param threshold the cut-off threshold for small elements in vector.
     */
    CSparseVector(const std::vector<double>& values, const std::vector<int32_t>& indexes, const int32_t nLength, const double threshold);

    /**
     Creates an empty sparse vector object with specific number of elements.

     @param nElements the number of elements in vector.
     @param threshold the cut-off thrshold for small matrix elements.
     */
    CSparseVector(const int32_t nElements, const double threshold);

    /**
     Creates a sparse vector object by copying other sparse vector object.

     @param source the sparse vector object.
     */
    CSparseVector(const CSparseVector& source);

    /**
     Creates a sparse vector object by moving other sparse vector object.

     @param source the sparse vector object.
     */
    CSparseVector(CSparseVector&& source) noexcept;

    /**
     Destroys a sparse vector object.
     */
    ~CSparseVector();

    /**
     Assigns a sparse vector object by copying other sparse vector object.

     @param source the sparse vector object.
     */
    CSparseVector& operator=(const CSparseVector& source);

    /**
     Assigns a sparse vector object by moving other sparse vector object.

     @param source the sparse vector object.
     */
    CSparseVector& operator=(CSparseVector&& source) noexcept;

    /**
     Compares sparse vector object with other sparse vector object.

     @param other the sparse vector object.
     @return true if sparse vector objects are equal, false otherwise.
     */
    bool operator==(const CSparseVector& other) const;

    /**
     Compares sparse vector object with other sparse vector object.

     @param other the sparse vector object.
     @return true if sparse vector objects are not equal, false otherwise.
     */
    bool operator!=(const CSparseVector& other) const;

    /**
     Gets constant pointer to sparse vector data.

     @return the constant pointer to first data element in vector matrix.
     */
    const double* values() const;

    /**
     Gets pointer to vector matrix data.

     @return the pointer to first data element in vector matrix.
     */
    double* values();

    /**
     Gets constant pointer to indexes of sparse vector.

     @return the constant pointer to first element of indexing vector.
     */
    const int32_t* indexes() const;

    /**
     Gets cut-off threshold of vector elements.

     @return the cut-off threshold.
     */
    double getThreshold() const;

    /**
     Computes sparsity number for sparse vector.

     @return the sparsity number.
     */
    double getSparsity() const;

    /**
     Converts sparse vector object to text output and insert it into output
     text stream.

     @param output the output text stream.
     @param source the sparse vector object.
     */
    friend std::ostream& operator<<(std::ostream& output, const CSparseVector& source);
};

#endif /* SparseVector_hpp */
