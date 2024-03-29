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

#ifndef ThreeIndexes_hpp
#define ThreeIndexes_hpp

#include <cstdint>
#include <ostream>

/**
 Class CThreeIndexes stores information about triple of indexes and provides
 functions to manipulate these indexes.

 @author Z. Rinkevicius
 */
class CThreeIndexes
{
    /**
     The first index from triple of indexes.
     */
    int32_t _iIndex;

    /**
     The second index from triple of indexes.
     */
    int32_t _jIndex;

    /**
     The third index from triple of indexes.
     */
    int32_t _kIndex;

   public:
    /**
     Creates an empty three indexes object.
     */
    CThreeIndexes();

    /**
     Creates a three indexes object.

     @param iIndex the first index from triple of indexes.
     @param jIndex the second index from triple of indexes.
     @param kIndex the third index from triple of indexes.
     */
    CThreeIndexes(const int32_t iIndex, const int32_t jIndex, const int32_t kIndex);

    /**
     Creates a three indexes object by copying other three indexes object.

     @param source the three indexes object.
     */
    CThreeIndexes(const CThreeIndexes& source);

    /**
     Creates a three indexes object by moving other three indexes object.

     @param source the three indexes object.
     */
    CThreeIndexes(CThreeIndexes&& source) noexcept;

    /**
     Destroys a three indexes object.
     */
    ~CThreeIndexes();

    /**
     Assigns a three indexes object by copying other three indexes object.

     @param source the three indexes object.
     */
    CThreeIndexes& operator=(const CThreeIndexes& source);

    /**
     Assigns a three indexes object by moving other three indexes object.

     @param source the three indexes object.
     */
    CThreeIndexes& operator=(CThreeIndexes&& source) noexcept;

    /**
     Compares three indexes object with other three indexes object.

     @param other the three indexes object.
     @return true if three indexes objects are equal, false otherwise.
     */
    bool operator==(const CThreeIndexes& other) const;

    /**
     Compares three indexes object with other three indexes object.

     @param other the three indexes object.
     @return true if three indexes objects are not equal, false otherwise.
     */
    bool operator!=(const CThreeIndexes& other) const;

    /**
     Gets first index from triple of indexes.

     @return the first index from triple of indexes.
     */
    int32_t first() const;

    /**
     Gets second index from triple of indexes.

     @return the second index from triple of indexes.
     */
    int32_t second() const;

    /**
     Gets third index from triple of indexes.

     @return the third index from triple of indexes.
     */
    int32_t third() const;

    /**
     Check if triple of indexes represents valid index triple (0..+Inf, +0..Inf,
     +0..Inf).

     @return true if triple of indexes is valid index triple, false - otherwise.
     */
    bool isValidTriple() const;

    /**
     Converts triple indexes object to text output and insert it into output text
     stream.

     @param output the output text stream.
     @param source the triple indexes object.
     */
    friend std::ostream& operator<<(std::ostream& output, const CThreeIndexes& source);
};

#endif /* ThreeIndexes_hpp */
