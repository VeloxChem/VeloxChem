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

#ifndef BasisFunction_hpp
#define BasisFunction_hpp

#include <cstddef>
#include <utility>
#include <vector>

/// @brief Class CBasisFunction stores data about single contracted GTO and
/// provides set of methods for handling of GTO data.
class CBasisFunction
{
   public:
    /// @brief The default constructor.
    CBasisFunction();

    /// @brief The constructor with exponents, normalization factors, and angular
    /// momentum.
    /// @param exponents The vector of exponents of primitive Gaussian functions.
    /// @param norms The vector of normalization factors of primitive Gaussian
    /// functions.
    /// @param angular_momentum The angular momentum of basis function.
    CBasisFunction(const std::vector<double> &exponents, const std::vector<double> &norms, const int angular_momentum);

    /// @brief The default copy constructor.
    /// @param other The basis function to be copied.
    CBasisFunction(const CBasisFunction &other);

    /// @brief The default move constructor.
    /// @param other The basis function to be moved.
    CBasisFunction(CBasisFunction &&other) noexcept;

    /// @brief The default destructor.
    ~CBasisFunction() = default;

    /// @brief The default copy assignment operator.
    /// @param other The basis function to be copy assigned.
    /// @return The assigned basis function.
    auto operator=(const CBasisFunction &other) -> CBasisFunction &;

    /// @brief The default move assignment operator.
    /// @param other The basis function to be move assigned.
    /// @return The assigned basis function.
    auto operator=(CBasisFunction &&other) noexcept -> CBasisFunction &;

    /// @brief The equality operator.
    /// @param other The basis function to be compared.
    /// @return True if basis functions are equal, False otherwise.
    auto operator==(const CBasisFunction &other) const -> bool;

    /// @brief The equality operator.
    /// @param other The basis function to be compared.
    /// @return True if basis functions are not equal, False otherwise.
    auto operator!=(const CBasisFunction &other) const -> bool;

    /// @brief Sets exponents of primittive Gaussian functions to specific vector
    /// of exponents.
    /// @param exponents The vector of exponents.
    auto set_exponents(const std::vector<double> &exponents) -> void;

    /// @brief Sets normalization factors of primitive Gaussian functions to
    /// specific vector of normalization factors.
    /// @param norms The vector of normalization factors.
    auto set_normalization_factors(const std::vector<double> &norms) -> void;

    /// @brief Set angular momentum of basis function.
    /// @param angular_momentum The angular momentum value to set.
    auto set_angular_momentum(const int angular_momentum) -> void;

    /// @brief Adds primittive Gaussian function to basis function.
    /// @param exponent The exponent of primitive Gaussian function.
    /// @param norm The normalization factor of primitive Gaussian function.
    auto add(const double exponent, const double norm) -> void;

    /// @brief Normalizes basis function.
    auto normalize() -> void;

    /// @brief Gets vector of exponents of primitive Gaussian functions.
    /// @return The vector of exponents of primitive Gaussian functions.
    auto get_exponents() const -> std::vector<double>;

    /// @brief Gets vector of normalization factors of primitive Gaussian
    /// functions.
    /// @return The vector of normalization factors of primitive Gaussian
    /// functions.
    auto get_normalization_factors() const -> std::vector<double>;

    /// @brief Gets angular momentum of basis function.
    /// @return The angular momentum of basis function.
    auto get_angular_momentum() const -> int;

    /// @brief Gets number of primitive Gaussian functions in basis function.
    /// @return The number of primitive Gaussian functions in basis function.
    auto number_of_primitive_functions() const -> size_t;

   private:
    /// @brief The vector of exponents of primitive Gaussian functions.
    std::vector<double> _exponents;

    /// @brief The vector of normalization factors of primitive Gaussian
    /// functions.
    std::vector<double> _norms;

    /// @brief The angular momentum of basis function.
    int _angular_momentum;

    /// @brief Rescales normalization factors to match normalization of spherical
    /// (l,0) component of basis function.
    auto _rescale() -> void;

    /// @brief Computes overlap between two primitive Gaussian functions.
    /// @param index The index pair of first and second primitve Gaussain
    /// functions.
    /// @return The overlap of two primitive Gaussian functions.
    auto _overlap(const std::pair<size_t, size_t> &index) const -> double;
};

#endif /* BasisFunction_hpp */
