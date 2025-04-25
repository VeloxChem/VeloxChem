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
