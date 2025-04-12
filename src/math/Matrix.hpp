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

#ifndef Matrix_hpp
#define Matrix_hpp

#include <cstdint>
#include <map>
#include <set>
#include <vector>

#include "MatrixType.hpp"
#include "SubMatrix.hpp"
#include "T2Pair.hpp"
#include "T4Index.hpp"

/**
 Class CMatrix stores matrix data and provides set of methods
 for handling of matrix data.

 @author Z. Rinkevicius
 */
class CMatrix
{
    /**
     The vector of submatrices.
     */
    std::map<T2Pair, CSubMatrix*> _sub_matrices;

    /**
     The type of matrix.
     */
    mat_t _mat_type;

    /**
     Gets set of row angular keys.
     */
    auto _getRowAngularKeys() const -> std::set<int64_t>;

    /**
     Gets set of column angular keys.
     */
    auto _getColumnAngularKeys() const -> std::set<int64_t>;

   public:
    /**
     Creates an empty matrix.
     */
    CMatrix();

    /**
     Creates a matrix.

     @param sub_matrices the map of submatrices.
     @param mat_type the matrix type.
     */
    CMatrix(const std::map<T2Pair, CSubMatrix>& sub_matrices, const mat_t mat_type);

    /**
     Creates a matrix.

     @param other the matrix to copy.
     */
    CMatrix(const CMatrix& other);

    /**
     Destroys a matrix.
     */
    ~CMatrix();

    /**
     Assigns a matrix  by copying other matrix.

     @param source the matrix to copy.
    */
    auto operator=(const CMatrix& source) -> CMatrix&;

    /**
     Adds submatrix to matrix.

     @param sub_matrix the submatrix to be added.
     @param angpair the angular pair of submatrix.
     */
    auto add(const CSubMatrix& sub_matrix, const T2Pair& angpair) -> void;

    /**
     Adds submatrix to matrix.

    @param dimensions the dimensions of submatrix to be added.
    @param angpair the angular pair of submatrix.
     */
    auto add(const T4Index& dimensions, const T2Pair& angpair) -> void;

    /**
     Set type  of matrix.

     @param mat_type the matrix type.
     */
    auto setType(const mat_t mat_type) -> void;

    /**
     Set matrix values to zero.
     */
    auto zero() -> void;

    /**
     Get vector of angular pairs from map  of submatrices.

     @return the vector of angular pairs.
     */
    auto getAngularPairs() const -> std::vector<T2Pair>;

    /**
     Get type  of matrix.

     @return the type of matrix.
     */
    auto getType() const -> mat_t;

    /**
     Get pointer to specific submatrix.

     @param angpair the angular pair of submatrix.
     @return the pointer to requested submatrix.
     */
    auto getSubMatrix(const T2Pair& angpair) -> CSubMatrix*;

    /**
     Get constant pointer to specific submatrix.

     @param angpair the angular pair of submatrix.
     @return the constant pointer to requested submatrix.
     */
    auto getSubMatrix(const T2Pair& angpair) const -> const CSubMatrix*;

    /**
     Checks if submatrix with requested angular pair is stored in matrix.

     @param angpair the angular pair of submatrix.
     @return True if the submatrix with this angular pair is stored, False otherwise.
     */
    auto isAngularOrder(const T2Pair& angpair) const -> bool;

    /**
     Gets number of rows in matrix.

     @return the number of rows in matrix.
     */
    auto getNumberOfRows() const -> int64_t;

    /**
     Gets number of columns in matrix.

     @return the number of columns in Fock matrix.
     */
    auto getNumberOfColumns() const -> int64_t;

    /**
     Reconstructs full matrix as submatrix with dimensions (0, 0, naos, naos).

     @return the full matrix.
     */
    auto getFullMatrix() const -> CSubMatrix;
};

#endif /* Matrix_hpp */
