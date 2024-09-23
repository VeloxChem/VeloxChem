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

#ifndef Matrices_hpp
#define Matrices_hpp

#include <map>
#include <string>
#include <vector>

#include "Matrix.hpp"

/// @brief Class CMatrices stores dictionary of matrices and provides set of methods
/// for handling its data.
class CMatrices
{
   public:
    /// @brief The default constructor.
    CMatrices();

    /// @brief The constructor with dictionary of matrices.
    /// @param matrices The map  of matrices.
    CMatrices(const std::map<std::string, CMatrix>& matrices);

    /// @brief The default copy constructor.
    /// @param other The dictionary of matrices  to be copied.
    CMatrices(const CMatrices& other);

    /// @brief The default move constructor.
    /// @param other The dictionary of matrices to be moved.
    CMatrices(CMatrices&& other) noexcept;

    /// @brief The default destructor.
    ~CMatrices();

    /// @brief The default copy assignment operator.
    /// @param other The dictionary of matrices  to be copy assigned.
    /// @return The assigned dictionary of matrices.
    auto operator=(const CMatrices& other) -> CMatrices&;

    /// @brief The default move assignment operator.
    /// @param other The dictionary of matrices  to be move assigned.
    /// @return The assigned dictionary of matrices.
    auto operator=(CMatrices&& other) noexcept -> CMatrices&;

    /// @brief The equality operator.
    /// @param other The dictionary of matrices to be compared.
    /// @return True if dictionaries of matrices are equal, False otherwise.
    auto operator==(const CMatrices& other) const -> bool;

    /// @brief The equality operator.
    /// @param other The dictionary of matrices to be compared.
    /// @return True if dictionaries of matrices are not equal, False otherwise.
    auto operator!=(const CMatrices& other) const -> bool;

    /// @brief Adds matrix to dictionary of matrices.
    /// @param matrix The matrix to be added.
    /// @param key The key of matrix to be added.
    auto add(const CMatrix& matrix, const std::string& key) -> void;

    /// @brief Adds matrix to dictionary of matrices.
    /// @param matrix The matrix to be added.
    /// @param key The key of matrix to be added.
    auto add(const CMatrix& matrix, const int key) -> void;

    /// @brief Set matrices values to zero in dictionary of matrices.
    auto zero() -> void;

    /// @brief Scales matrices values by factor in dictionary of matrices.
    /// @param factor The factor to scale matrices values.
    auto scale(const double factor) -> void;

    /// @brief Symmetrizes matrices in dictionary of matrices.
    auto symmetrize() -> void;

    /// @brief Get vector of keys from dictionary of matrices.
    /// @return The vector of keys.
    auto keys() const -> std::vector<std::string>;

    /// @brief Gets pointer to requested matrix.
    /// @param key The key of requested matrix.
    /// @return The pointer to matrix.
    auto matrix(const std::string& key) -> CMatrix*;

    /// @brief Gets pointer to requested matrix.
    /// @param key The key of requested matrix.
    /// @return The pointer to matrix.
    auto matrix(const std::string& key) const -> const CMatrix*;

    /// @brief Gets constant pointer to requested matrix.
    /// @param key The key of requested matrix.
    /// @return The constant pointer to matrix.
    auto matrix(const int key) -> CMatrix*;

    /// @brief Gets constant pointer to requested matrix.
    /// @param key The key of requested matrix.
    /// @return The constant pointer to matrix.
    auto matrix(const int key) const -> const CMatrix*;

   private:
    /// @brief The map of matrices.
    std::map<std::string, CMatrix*> _matrices;

    /// @brief Deallocates matrices in dictionary of matrices.
    auto _deallocate() -> void;
};

#endif /* Matrices_hpp */
