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

#ifndef GtoBlock_hpp
#define GtoBlock_hpp

#include <vector>

#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "Point.hpp"

/// @brief Class CGtoBlock stores data about basis functions block and provides
/// set of methods for manipulating with basis functions block.
class CGtoBlock
{
   public:
    /// @brief The default constructor.
    CGtoBlock();

    /// @brief The constructor with basis functions data.
    /// @param coordinates The vector of basis function Cartesian coordinates.
    /// @param exponents The vector of primitive basis function exponents.
    /// @param norms The vector of primitive basis function normalization
    /// factors.
    /// @param orb_indices The vector of atomic orbital indices.
    /// @param atm_indices The vector of atomic indices.
    /// @param angular_momentum The angular momentum of basis functions.
    /// @param npgtos The number of primitive basis functions in basis function.
    CGtoBlock(const std::vector<TPoint<double>> &coordinates,
              const std::vector<double>         &exponents,
              const std::vector<double>         &norms,
              const std::vector<size_t>         &orb_indices,
              const std::vector<int>            &atm_indices,
              const int                          angular_momentum,
              const int                          npgtos);

    /// @brief The constructor molecular basis and molecule.
    /// @param basis The molecular basis.
    /// @param molecule The molecule.
    /// @param angular_momentum The angular momentum of basis functions.
    /// @param npgtos The number of primitive basis functions in basis function.
    CGtoBlock(const CMolecularBasis &basis, const CMolecule &molecule, const int angular_momentum, const int npgtos);

    /// @brief The constructor molecular basis and molecule for selected atoms in
    /// molecule.
    /// @param basis The molecular basis.
    /// @param molecule The molecule.
    /// @param atoms The vector of atoms to select.
    /// @param angular_momentum The angular momentum of basis functions.
    /// @param npgtos The number of primitive basis functions in basis function.
    CGtoBlock(const CMolecularBasis &basis, const CMolecule &molecule, const std::vector<int> &atoms, const int angular_momentum, const int npgtos);

    /// @brief The default copy constructor.
    /// @param other The basis functions block to be copied.
    CGtoBlock(const CGtoBlock &other);

    /// @brief The default move constructor.
    /// @param other The basis functions block to be moved.
    CGtoBlock(CGtoBlock &&other) noexcept;

    /// @brief The default destructor.
    ~CGtoBlock() = default;

    /// @brief The default copy assignment operator.
    /// @param other The basis functions block to be copy assigned.
    /// @return The assigned basis functions block.
    auto operator=(const CGtoBlock &other) -> CGtoBlock &;

    /// @brief The default move assignment operator.
    /// @param other The basis functions block to be move assigned.
    /// @return The assigned basis functions block.
    auto operator=(CGtoBlock &&other) noexcept -> CGtoBlock &;

    /// @brief The equality operator.
    /// @param other The basis functions block to be compared.
    /// @return True if basis functions blocks are equal, False otherwise.
    auto operator==(const CGtoBlock &other) const -> bool;

    /// @brief The equality operator.
    /// @param other The basis functions block to be compared.
    /// @return True if basis functions blocks are not equal, False otherwise.
    auto operator!=(const CGtoBlock &other) const -> bool;

    /// @brief Gets vector of basis function Cartesian coordinates.
    /// @return The vector of Cartesian coordinates.
    auto coordinates() const -> std::vector<TPoint<double>>;
    
    /// @brief Gets vector of basis function Cartesian coordinates.
    /// @param origin The origin of Cartesian coordinates system.
    /// @return The vector of Cartesian coordinates.
    auto coordinates(const TPoint<double>& origin) const -> std::vector<TPoint<double>>;

    /// @brief Gets vector of basis functions exponents.
    /// @return The vector of exponents.
    auto exponents() const -> std::vector<double>;

    /// @brief Gets vector of basis functions normalization factors.
    /// @return The vector of normalization factors.
    auto normalization_factors() const -> std::vector<double>;

    /// @brief Gets vector of orbital indices of basis functions.
    /// @return The vector of orbital indices.
    auto orbital_indices() const -> std::vector<size_t>;

    /// @brief Gets vector of atomic indices of basis functions.
    /// @return The vector of atomic indices.
    auto atomic_indices() const -> std::vector<int>;

    /**
     Gets vector of atomic orbitals indexes of contracted GTOs.

     @return the vector of atomic orbitals indexes of GTOs.
     */
    auto getAtomicOrbitalsIndexes() const -> std::vector<int>;

    auto getAtomicOrbitalsIndexesForCartesian(const int ncgtos_d = 0) const -> std::vector<int>;

    auto getCartesianToSphericalMappingForP() const -> std::unordered_map<int, std::vector<std::pair<int, double>>>;

    auto getCartesianToSphericalMappingForD() const -> std::unordered_map<int, std::vector<std::pair<int, double>>>;

    auto getCartesianToSphericalMappingForF(const int ncgtos_d) const -> std::unordered_map<int, std::vector<std::pair<int, double>>>;

    /// @brief Gets angular momentum of basis functions.
    /// @return The angular momentum of basis functionss.
    auto angular_momentum() const -> int;

    /// @brief Gets number of primitive basis functions in basis function.
    /// @return The number of primitive basis functions in basis function.
    auto number_of_primitives() const -> int;

    /// @brief Gets number of basis functions in basis functions block.
    /// @return The number of basis functions in basis functions block.
    auto number_of_basis_functions() const -> int;

   private:
    /// @brief The vector of Cartesian coordinates of basis functions.
    std::vector<TPoint<double>> _coordinates;

    /// @brief The vector of exponents of primitive basis functions.
    std::vector<double> _exponents;

    /// @brief The vector of normalization factors of primitive basis functions.
    std::vector<double> _norms;

    /// @brief The vector of atomic orbitals indices of basis functions.
    std::vector<size_t> _orb_indices;

    /// @brief The vector of atomic indices of basis functions.
    std::vector<int> _atm_indices;

    /// @brief The angular momentum of basis functions.
    int _angular_momentum;

    /// @brief The number of primitive basis functions in basis function.
    int _npgtos;
};

#endif /* GtoBlock_hpp */
