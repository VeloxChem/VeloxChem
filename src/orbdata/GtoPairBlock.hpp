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

#ifndef GtoPairBlock_hpp
#define GtoPairBlock_hpp

#include <array>
#include <string>
#include <utility>
#include <vector>

#include "GtoBlock.hpp"
#include "Point.hpp"

/// @brief Class CGtoPairBlock stores data about basis function pairs block and
/// provides set of methods for manipulating with basis function pairs block.
class CGtoPairBlock
{
   public:
    /// Creates an empty basis function pairs block.
    CGtoPairBlock() = default;

    /// @brief Creates a basis function pairs block.
    /// @param bra_coordinates The vector of basis functions Cartesian
    /// coordinates on bra side.
    /// @param ket_coordinates The vector of basis functions Cartesian
    /// coordinates on ket side.
    /// @param bra_exponents The vector of exponents of primitive basis functions
    /// on bra side.
    /// @param ket_exponents The vector of exponents of primitive basis functions
    /// on ket side.
    /// @param norms The vector of normalization factors of primitive basis
    /// function pairs.
    /// @param overlaps The vector of overlap factors of primitive basis function
    /// pairs.
    /// @param bra_orb_indices The vector of  AO indices on bra side.
    /// @param ket_orb_indices The vector of  AO indices on ket side.
    /// @param bra_atm_indices The vector of  atomic indices on bra side.
    /// @param ket_atm_indices The vector of  atomic indices on ket side.
    /// @param angular_momentums The angular momentums of basis function pair.
    /// @param nppairs The number of primitive basis function pairs in contracted
    /// GTO pairs.
    CGtoPairBlock(const std::vector<TPoint<double>> &bra_coordinates,
                  const std::vector<TPoint<double>> &ket_coordinates,
                  const std::vector<double>         &bra_exponents,
                  const std::vector<double>         &ket_exponents,
                  const std::vector<double>         &norms,
                  const std::vector<double>         &overlaps,
                  const std::vector<size_t>         &bra_orb_indices,
                  const std::vector<size_t>         &ket_orb_indices,
                  const std::vector<int>            &bra_atm_indices,
                  const std::vector<int>            &ket_atm_indices,
                  const std::pair<int, int>         &angular_momentums,
                  const int                          nppairs);

    /// @brief Creates a basis function pairs block.
    /// @param gto_block The basis functions block on bra and ket side.
    CGtoPairBlock(const CGtoBlock &gto_block);

    /// @brief Creates a basis function pairs block.
    /// @param bra_gto_block The basis functions block on bra side.
    /// @param ket_gto_block The basis function block on ket side.
    CGtoPairBlock(const CGtoBlock &bra_gto_block, const CGtoBlock &ket_gto_block);

    /// @brief The default copy constructor.
    /// @param other The basis function pairs block to be copied.
    CGtoPairBlock(const CGtoPairBlock &other);

    /// @brief The default move constructor.
    /// @param other The basis functions pairs block to be moved.
    CGtoPairBlock(CGtoPairBlock &&other) noexcept;

    /// @brief The default destructor.
    ~CGtoPairBlock() = default;

    /// @brief The default copy assignment operator.
    /// @param other The basis function pairs block to be copy assigned.
    /// @return The assigned basis function pairs block.
    auto operator=(const CGtoPairBlock &other) -> CGtoPairBlock &;

    /// @brief The default move assignment operator.
    /// @param other The basis function pairs block to be move assigned.
    /// @return The assigned basis function pairs block.
    auto operator=(CGtoPairBlock &&other) noexcept -> CGtoPairBlock &;

    /// @brief The equality operator.
    /// @param other The basis function pairs block to be compared.
    /// @return True if basis function pairs blocks are equal, False otherwise.
    auto operator==(const CGtoPairBlock &other) const -> bool;

    /// @brief The equality operator.
    /// @param other The basis function pairs  block to be compared.
    /// @return True if basis function pairs blocks are not equal, False
    /// otherwise.
    auto operator!=(const CGtoPairBlock &other) const -> bool;

    /// @brief Gets vector of GTO pair Cartesian coordinates on bra side.
    /// @return The vector of Cartesian coordinates.
    auto bra_coordinates() const -> std::vector<TPoint<double>>;

    /// @brief Gets vector of GTO pair Cartesian coordinates on ket side.
    /// @return The vector of Cartesian coordinates.
    auto ket_coordinates() const -> std::vector<TPoint<double>>;

    /// @brief Gets vector of primitive basis function exponents on bra side.
    /// @return The vector of primtive basis function exponentns.
    auto bra_exponents() const -> std::vector<double>;

    /// @brief Gets vector of primitive basis function exponents on ket side.
    /// @return The vector of primtive basis function exponentns.
    auto ket_exponents() const -> std::vector<double>;

    /// @brief Gets vector of primitive basis function pairs normalization
    /// factors.
    /// @return The vector of normalization factors.
    auto normalization_factors() const -> std::vector<double>;

    /// @brief Gets vector of primitive basis function pairs overlap factors.
    /// @return The vector of overlap factors.
    auto overlap_factors() const -> std::vector<double>;

    /// @brief Gets vector of orbital indices of basis function pairs on bra
    /// side.
    /// @return The vector of orbital indices.
    auto bra_orbital_indices() const -> std::vector<size_t>;

    /// @brief Gets vector of orbital indices of basis function pairs on ket
    /// side.
    /// @return The vector of orbital indices.
    auto ket_orbital_indices() const -> std::vector<size_t>;

    /// @brief Gets vector of atomic indices of basis function pairs on bra side.
    /// @return The vector of atomic indices.
    auto bra_atomic_indices() const -> std::vector<int>;

    /// @brief Gets vector of atomic indices of basis function pairs on ket side.
    /// @return The vector of atomic indices.
    auto ket_atomic_indices() const -> std::vector<int>;

    /// @brief Gets angular momentums of basis function pair.
    /// @return The angular momentum of basis function pair.
    auto angular_momentums() const -> std::pair<int, int>;

    /// @brief Gets number of primitive basis function pairs in basis function
    /// pair.
    /// @return The number of primtive basis function pairs.
    auto number_of_primitive_pairs() const -> int;

    /// @brief Gets number of basis function pairs in basis function pairs block.
    /// @return The number of basis function pairs.
    auto number_of_contracted_pairs() const -> size_t;
    
    /// @brief Gets number of unique terms (integral bra or ket side of integral) generated by this basis function pairs block.
    /// @return The number of unique terms.
    auto unique_terms() const -> size_t;

   private:
    /// @brief The vector of Cartesian coordinates of contracted GTO pairs on bra
    /// side.
    std::vector<TPoint<double>> _bra_coordinates;

    /// @brief The vector of Cartesian coordinates of contracted GTO pairs on ket
    /// side.
    std::vector<TPoint<double>> _ket_coordinates;

    /// @brief The vector of exponents of primitive GTOs pairs on bra side.
    std::vector<double> _bra_exponents;

    /// @brief The vector of exponents of primitive GTOs pairs on ket side.
    std::vector<double> _ket_exponents;

    /// @brief The vector of normalization factors of primitive GTO pairs.
    std::vector<double> _norms;

    /// @brief The vector of overlap factors of primitive GTO pairs.
    std::vector<double> _overlaps;

    /// @brief The vector of AO indices of contracted GTO pairs on bra side.
    std::vector<size_t> _bra_orb_indices;

    /// @brief The vector of AO indices of contracted GTO pairs on ket side.
    std::vector<size_t> _ket_orb_indices;

    /// @brief The vector of atomic indices of contracted GTO pairs on bra side.
    std::vector<int> _bra_atm_indices;

    /// @brief The vector of atomic indices of contracted GTO pairs on ket side.
    std::vector<int> _ket_atm_indices;

    /// @brief The angular momentums of contracted GTO pair.
    std::pair<int, int> _angular_momentums;

    /// @brief The number of primitive GTOs in contracted GTO pair.
    int _nppairs;
};

#endif /* GtoPairBlock_hpp */
