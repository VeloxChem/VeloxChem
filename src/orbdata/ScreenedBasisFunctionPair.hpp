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

#ifndef ScreenedBasisFunctionPair_hpp
#define ScreenedBasisFunctionPair_hpp

#include <cstddef>
#include <utility>
#include <vector>

#include "BasisFunction.hpp"
#include "Point.hpp"

/// @brief Class CScreenedBasisFunctionPair stores data about a screened pair of
/// basis functions: the two basis functions together with the Cartesian
/// coordinates and atomic orbital indices of the bra and ket atom pairs that
/// survive a screening predicate. Coordinates are stored as separate x, y, z
/// vectors (structure of arrays) on the bra and ket sides.
class CScreenedBasisFunctionPair
{
   public:
    /// @brief The default constructor.
    CScreenedBasisFunctionPair();

    /// @brief The constructor with explicit screened pair data.
    /// @param bra_function The basis function on bra side.
    /// @param ket_function The basis function on ket side.
    /// @param bra_x The vector of bra atom Cartesian X coordinates.
    /// @param bra_y The vector of bra atom Cartesian Y coordinates.
    /// @param bra_z The vector of bra atom Cartesian Z coordinates.
    /// @param ket_x The vector of ket atom Cartesian X coordinates.
    /// @param ket_y The vector of ket atom Cartesian Y coordinates.
    /// @param ket_z The vector of ket atom Cartesian Z coordinates.
    /// @param orb_indices The vector of (bra, ket) atomic orbital index pairs.
    CScreenedBasisFunctionPair(const CBasisFunction                         &bra_function,
                               const CBasisFunction                         &ket_function,
                               const std::vector<double>                    &bra_x,
                               const std::vector<double>                    &bra_y,
                               const std::vector<double>                    &bra_z,
                               const std::vector<double>                    &ket_x,
                               const std::vector<double>                    &ket_y,
                               const std::vector<double>                    &ket_z,
                               const std::vector<std::pair<size_t, size_t>> &orb_indices);

    /// @brief The constructor screening symmetric atom pairs of a single basis
    /// function. The ket basis function is identical to the bra basis function
    /// and the atom pairs are generated in a triangular pattern (i <= j),
    /// including the diagonal self-pairs.
    /// @param function The basis function on bra and ket side.
    /// @param coordinates The vector of atom Cartesian coordinates.
    /// @param orb_indices The vector of atomic orbital indices.
    /// @param keep_pair The predicate selecting important atom pairs; called as
    /// keep_pair(bra_center, ket_center) and returning true to keep the pair.
    template <typename Predicate>
    CScreenedBasisFunctionPair(const CBasisFunction              &function,
                               const std::vector<TPoint<double>> &coordinates,
                               const std::vector<size_t>         &orb_indices,
                               Predicate                        &&keep_pair)

        : _bra_function(function)

        , _ket_function(function)

        , _bra_x{}
        , _bra_y{}
        , _bra_z{}
        , _ket_x{}
        , _ket_y{}
        , _ket_z{}
        , _orb_indices{}
    {
        const auto natoms = coordinates.size();

        for (size_t i = 0; i < natoms; i++)
        {
            for (size_t j = i; j < natoms; j++)
            {
                if (keep_pair(coordinates[i], coordinates[j]))
                {
                    _append(coordinates[i], coordinates[j], orb_indices[i], orb_indices[j]);
                }
            }
        }
    }

    /// @brief The constructor screening atom pairs of two distinct basis
    /// functions. The atom pairs are generated as the full product of bra and
    /// ket atoms.
    /// @param bra_function The basis function on bra side.
    /// @param ket_function The basis function on ket side.
    /// @param bra_coordinates The vector of bra atom Cartesian coordinates.
    /// @param bra_orb_indices The vector of bra atomic orbital indices.
    /// @param ket_coordinates The vector of ket atom Cartesian coordinates.
    /// @param ket_orb_indices The vector of ket atomic orbital indices.
    /// @param keep_pair The predicate selecting important atom pairs; called as
    /// keep_pair(bra_center, ket_center) and returning true to keep the pair.
    template <typename Predicate>
    CScreenedBasisFunctionPair(const CBasisFunction              &bra_function,
                               const CBasisFunction              &ket_function,
                               const std::vector<TPoint<double>> &bra_coordinates,
                               const std::vector<size_t>         &bra_orb_indices,
                               const std::vector<TPoint<double>> &ket_coordinates,
                               const std::vector<size_t>         &ket_orb_indices,
                               Predicate                        &&keep_pair)

        : _bra_function(bra_function)

        , _ket_function(ket_function)

        , _bra_x{}
        , _bra_y{}
        , _bra_z{}
        , _ket_x{}
        , _ket_y{}
        , _ket_z{}
        , _orb_indices{}
    {
        const auto nbra = bra_coordinates.size();

        const auto nket = ket_coordinates.size();

        for (size_t i = 0; i < nbra; i++)
        {
            for (size_t j = 0; j < nket; j++)
            {
                if (keep_pair(bra_coordinates[i], ket_coordinates[j]))
                {
                    _append(bra_coordinates[i], ket_coordinates[j], bra_orb_indices[i], ket_orb_indices[j]);
                }
            }
        }
    }

    /// @brief The default copy constructor.
    /// @param other The screened basis function pair to be copied.
    CScreenedBasisFunctionPair(const CScreenedBasisFunctionPair &other);

    /// @brief The default move constructor.
    /// @param other The screened basis function pair to be moved.
    CScreenedBasisFunctionPair(CScreenedBasisFunctionPair &&other) noexcept;

    /// @brief The default destructor.
    ~CScreenedBasisFunctionPair() = default;

    /// @brief The default copy assignment operator.
    /// @param other The screened basis function pair to be copy assigned.
    /// @return The assigned screened basis function pair.
    auto operator=(const CScreenedBasisFunctionPair &other) -> CScreenedBasisFunctionPair &;

    /// @brief The default move assignment operator.
    /// @param other The screened basis function pair to be move assigned.
    /// @return The assigned screened basis function pair.
    auto operator=(CScreenedBasisFunctionPair &&other) noexcept -> CScreenedBasisFunctionPair &;

    /// @brief The equality operator.
    /// @param other The screened basis function pair to be compared.
    /// @return True if screened basis function pairs are equal, False otherwise.
    auto operator==(const CScreenedBasisFunctionPair &other) const -> bool;

    /// @brief The non-equality operator.
    /// @param other The screened basis function pair to be compared.
    /// @return True if screened basis function pairs are not equal, False
    /// otherwise.
    auto operator!=(const CScreenedBasisFunctionPair &other) const -> bool;

    /// @brief Gets the basis function on bra side.
    /// @return The bra basis function.
    auto bra_function() const -> CBasisFunction;

    /// @brief Gets the basis function on ket side.
    /// @return The ket basis function.
    auto ket_function() const -> CBasisFunction;

    /// @brief Gets the vector of bra atom Cartesian X coordinates.
    /// @return The vector of X coordinates.
    auto bra_x() const -> std::vector<double>;

    /// @brief Gets the vector of bra atom Cartesian Y coordinates.
    /// @return The vector of Y coordinates.
    auto bra_y() const -> std::vector<double>;

    /// @brief Gets the vector of bra atom Cartesian Z coordinates.
    /// @return The vector of Z coordinates.
    auto bra_z() const -> std::vector<double>;

    /// @brief Gets the vector of ket atom Cartesian X coordinates.
    /// @return The vector of X coordinates.
    auto ket_x() const -> std::vector<double>;

    /// @brief Gets the vector of ket atom Cartesian Y coordinates.
    /// @return The vector of Y coordinates.
    auto ket_y() const -> std::vector<double>;

    /// @brief Gets the vector of ket atom Cartesian Z coordinates.
    /// @return The vector of Z coordinates.
    auto ket_z() const -> std::vector<double>;

    /// @brief Gets the vector of (bra, ket) atomic orbital index pairs.
    /// @return The vector of orbital index pairs.
    auto orbital_indices() const -> std::vector<std::pair<size_t, size_t>>;

    /// @brief Gets the number of screened basis function pairs.
    /// @return The number of pairs.
    auto number_of_pairs() const -> size_t;

   private:
    /// @brief Appends a single screened atom pair to the structure of arrays.
    /// @param bra_center The Cartesian coordinates of the bra atom.
    /// @param ket_center The Cartesian coordinates of the ket atom.
    /// @param bra_index The bra atomic orbital index.
    /// @param ket_index The ket atomic orbital index.
    auto _append(const TPoint<double> &bra_center, const TPoint<double> &ket_center, const size_t bra_index, const size_t ket_index) -> void;

    /// @brief The basis function on bra side.
    CBasisFunction _bra_function;

    /// @brief The basis function on ket side.
    CBasisFunction _ket_function;

    /// @brief The vector of bra atom Cartesian X coordinates.
    std::vector<double> _bra_x;

    /// @brief The vector of bra atom Cartesian Y coordinates.
    std::vector<double> _bra_y;

    /// @brief The vector of bra atom Cartesian Z coordinates.
    std::vector<double> _bra_z;

    /// @brief The vector of ket atom Cartesian X coordinates.
    std::vector<double> _ket_x;

    /// @brief The vector of ket atom Cartesian Y coordinates.
    std::vector<double> _ket_y;

    /// @brief The vector of ket atom Cartesian Z coordinates.
    std::vector<double> _ket_z;

    /// @brief The vector of (bra, ket) atomic orbital index pairs.
    std::vector<std::pair<size_t, size_t>> _orb_indices;
};

#endif /* ScreenedBasisFunctionPair_hpp */
