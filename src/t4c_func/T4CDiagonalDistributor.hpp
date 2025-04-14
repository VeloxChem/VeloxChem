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

#ifndef T4CDiagonalDistributor_hpp
#define T4CDiagonalDistributor_hpp

#include <cstddef>
#include <utility>
#include <vector>

#include "SimdArray.hpp"

/// @brief Class CT4CDiagonalDistributor provides methods for distributing vector of maximum values of integral shells.
class CT4CDiagonalDistributor
{
   public:
    /// @brief Creates a maximum values of integral shells distributor.
    /// @param max_values  The pointer to vector of maximum values.
    CT4CDiagonalDistributor(double* max_values);

    /// @brief The default copy constructor.
    /// @param other The distributor to be copied.
    CT4CDiagonalDistributor(const CT4CDiagonalDistributor& other) = delete;

    /// @brief The default move constructor.
    /// @param other The distributor to be moved.
    CT4CDiagonalDistributor(CT4CDiagonalDistributor&& other) noexcept = delete;

    /// @brief The default destructor.
    ~CT4CDiagonalDistributor() = default;

    /// @brief The default copy assignment operator.
    /// @param other The distributor to be copy assigned.
    /// @return The assigned distributor.
    auto operator=(const CT4CDiagonalDistributor& other) -> CT4CDiagonalDistributor& = delete;

    /// @brief The default move assignment operator.
    /// @param other The distributor to be move assigned.
    /// @return The assigned distributor.
    auto operator=(CT4CDiagonalDistributor&& other) noexcept -> CT4CDiagonalDistributor& = delete;

    /// @brief The equality operator.
    /// @param other The distributor to be compared.
    /// @return True if distributors are equal, False otherwise.
    auto operator==(const CT4CDiagonalDistributor& other) const -> bool = delete;

    /// @brief The equality operator.
    /// @param other The distributor to be compared.
    /// @return True if distributors are not equal, False otherwise.
    auto operator!=(const CT4CDiagonalDistributor& other) const -> bool = delete;

    /// @brief Distributes maximum values of integral shells  into vector.
    /// @param max_integrals The vector of integral values.
    /// @param gto_range The index of the range [ket_first, ket_last) of integrals.
    auto distribute(const std::vector<double>& max_integrals, const std::pair<size_t, size_t>& gto_range) -> void;

   private:
    /// @brief The vector of maximum values of integral shells.
    double* _max_values;
};

#endif /* T4CDiagonalDistributor_hpp */
