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

#ifndef Dense4DTensor_hpp
#define Dense4DTensor_hpp

#include <vector>

/**
 Class CDense4DTensor stores 4D tensor in coordinate format (zero-based
 indexing scheme).
 */
class CDense4DTensor
{
    /**
     The dimension of the first index.
     */
    int _iIndex;

    /**
     The dimension of the second index.
     */
    int _jIndex;

    /**
     The dimension of the third index.
     */
    int _kIndex;

    /**
     The dimension of the fourth index.
     */
    int _lIndex;

    /**
     The vector of tensor element values.
     */
    std::vector<double> _values;

   public:
    /**
     Creates an empty dense matrix object.
     */
    CDense4DTensor();

    /**
     Creates a dense matrix object from vector of values.

     @param values the vector of tensor elements.
     @param iIndex the dimension of the first index.
     @param jIndex the dimension of the second index.
     @param kIndex the dimension of the third index.
     @param lIndex the dimension of the fourth index.
     */
    CDense4DTensor(const std::vector<double>& values, const int iIndex, const int jIndex, const int kIndex, const int lIndex);

    /**
     Creates an empty tensor object with specific dimensions.

     @param iIndex the dimension of the first index.
     @param jIndex the dimension of the second index.
     @param kIndex the dimension of the third index.
     @param lIndex the dimension of the fourth index.
     */
    CDense4DTensor(const int iIndex, const int jIndex, const int kIndex, const int lIndex);

    /**
     Creates an empty tensor object with all dimensions identical.

     @param nRows the number of rows in matrix.
     */
    CDense4DTensor(const int nRows);

    /**
     Creates a 4D tensor object by copying other 4D tensor object.

     @param source the 4D tensor object.
     */
    CDense4DTensor(const CDense4DTensor& source);

    /**
     Creates a 4D tensor object by moving other 4D tensor object.

     @param source the 4D tensor object.
     */
    CDense4DTensor(CDense4DTensor&& source) noexcept;

    /**
     Destroys a dense matrix object.
     */
    ~CDense4DTensor();

    /**
     Assigns 4D tensor object by copying other 4D tensor object.

     @param source the 4D tensor object.
     */
    CDense4DTensor& operator=(const CDense4DTensor& source);

    /**
     Assigns a 4D tensor object by moving other 4D tensor object.

     @param source the 4D tensor object.
     */
    CDense4DTensor& operator=(CDense4DTensor&& source) noexcept;

    /**
     Compares 4D tensor object with other 4D tensor object.

     @param other the 4D tensor object.
     @return true if 4D tensor objects are equal, false otherwise.
     */
    bool operator==(const CDense4DTensor& other) const;

    /**
     Compares 4D tensor object with other 4D tensor object.

     @param other the other 4D tensor object.
     @return true if 4D tensor objects are not equal, false otherwise.
     */
    bool operator!=(const CDense4DTensor& other) const;

    /**
     Sets all values in 4D tensor to zero.
     */
    void zero();

    /**
     Gets the dimension of the first index of the 4D tensor.

     @return the first index.
     */
    int getiIndex() const;

    /**
     Gets the dimension of the second index of the 4D tensor.

     @return the second index.
     */
    int getjIndex() const;

    /**
     Gets the dimension of the third index of the 4D tensor.

     @return the third index.
     */
    int getkIndex() const;

    /**
     Gets the dimension of the forth index of the 4D tensor.

     @return the fourth index.
     */
    int getlIndex() const;

    /**
     Gets number of elements in 4D tensor.

     @return the number of elements.
     */
    int getNumberOfElements() const;

    /**
     Gets constant pointer to first element of 4D tensor.

     @return the constant pointer to first element of dense matrix.
     */
    const double* values() const;

    /**
     Gets pointer to first element of 4D tensor.

     @return the pointer to first element of 4D tensor.
     */
    double* values();
};

#endif /* Dense4DTensor_hpp */
