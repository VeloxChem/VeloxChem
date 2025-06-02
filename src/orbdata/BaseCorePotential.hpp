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

#ifndef BaseCorePotential_hpp
#define BaseCorePotential_hpp

#include <cstddef>
#include <vector>

/// @brief Class CBaseCorePotential stores data about base core potential and
/// provides set of methods for handling of base core potential data.
class CBaseCorePotential
{
   public:
    /// @brief The default constructor.
    CBaseCorePotential();

    /// @brief The constructor with exponents, expansion factors, and radial orders.
    /// @param exponents The vector of exponents of primitive local potentials.
    /// @param factors The vector of expansion factors of primitive local potential.
    /// @param radial_orders The vector of radial orders of primitive local potentials.
    CBaseCorePotential(const std::vector<double> &exponents,
                       const std::vector<double> &factors,
                       const std::vector<int>    &radial_orders);

    /// @brief The default copy constructor.
    /// @param other The base core potential to be copied.
    CBaseCorePotential(const CBaseCorePotential &other);

    /// @brief The default move constructor.
    /// @param other The base core potential to be moved.
    CBaseCorePotential(CBaseCorePotential &&other) noexcept;

    /// @brief The default destructor.
    ~CBaseCorePotential() = default;

    /// @brief The default copy assignment operator.
    /// @param other The base core potential to be copy assigned.
    /// @return The assigned base core potential.
    auto operator=(const CBaseCorePotential &other) -> CBaseCorePotential&;

    /// @brief The default move assignment operator.
    /// @param other The base core potential to be move assigned.
    /// @return The assigned base core potential.
    auto operator=(CBaseCorePotential &&other) noexcept -> CBaseCorePotential &;

    /// @brief The equality operator.
    /// @param other The base core potential to be compared.
    /// @return True if base core potentials are equal, False otherwise.
    auto operator==(const CBaseCorePotential &other) const -> bool;

    /// @brief The equality operator.
    /// @param other The base core potential to be compared.
    /// @return True if base core potentials are not equal, False otherwise.
    auto operator!=(const CBaseCorePotential &other) const -> bool;

    /// @brief Sets exponents of primittive base core potentials to specific vector
    /// of exponents.
    /// @param exponents The vector of exponents.
    auto set_exponents(const std::vector<double> &exponents) -> void;

    /// @brief Sets expansion factors of primitive base core potentials to
    /// specific vector of expansion factors.
    /// @param factors The vector of expansion factors.
    auto set_factors(const std::vector<double> &factors) -> void;
    
    /// @brief Sets radial orders of primitive base core potentials to
    /// specific vector of radial orders.
    /// @param radial_orders The vector of radial orders.
    auto set_radial_orders(const std::vector<int> &radial_orders) -> void;

    /// @brief Adds primittive base core potential to base core potential.
    /// @param exponent The exponent of primitive base core potential.
    /// @param factor The expansion factor of primitive base core potential.
    /// @param radial_order  The radial order of primitive base core potential.
    auto add(const double exponent, const double factor, const int radial_order) -> void;

    /// @brief Gets vector of exponents of primitive base core potentials.
    /// @return The vector of exponents of primitive base core potentials.
    auto get_exponents() const -> std::vector<double>;

    /// @brief Gets vector of expansion factors of primitive base core potentials.
    /// @return The vector of expansion factors of primitive base core potentials.
    auto get_factors() const -> std::vector<double>;

    /// @brief Gets vector of radial orders of primitive base core potentials.
    /// @return The vector of radial orders of primitive base core potentials.
    auto get_radial_orders() const -> std::vector<int>;

    /// @brief Gets number of primitive  base core potentials in base core potential.
    /// @return The number of primitive  base core potentials in base core potential.
    auto number_of_primitive_potentials() const -> size_t;
    
    /// @brief Checks if radial orders are compatable with analytical ECP integration.
    /// @return True if radial orders compatable with analytical ECP integration, false otherwise.
    auto is_valid_radial_orders() const -> bool;

   private:
    /// @brief The vector of exponents of primitive local potentials.
    std::vector<double> _exponents;

    /// @brief The vector of expansion factors of primitive local potentials.
    std::vector<double> _factors;

    /// @brief The vector of radial orders of primitive local potentials.
    std::vector<int> _radial_orders;
};

#endif /* BaseCorePotential_hpp */
