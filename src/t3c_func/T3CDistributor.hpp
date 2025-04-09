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

#ifndef T3CDistributor_hpp
#define T3CDistributor_hpp

#include <cstddef>
#include <utility>
#include <vector>

#include "SimdArray.hpp"
#include "T3FlatBuffer.hpp"

/// @brief Class CT3CDistributor provides methods for distributing vector of into flat buffer.
class CT3CDistributor
{
   public:
    /// @brief Creates an integral shells distributor.
    /// @param values  The pointer to flat tensor.
    CT3CDistributor(CT3FlatBuffer<double>* values);

    /// @brief The default copy constructor.
    /// @param other The distributor to be copied.
    CT3CDistributor(const CT3CDistributor& other) = delete;

    /// @brief The default move constructor.
    /// @param other The distributor to be moved.
    CT3CDistributor(CT3CDistributor&& other) noexcept = delete;

    /// @brief The default destructor.
    ~CT3CDistributor() = default;

    /// @brief The default copy assignment operator.
    /// @param other The distributor to be copy assigned.
    /// @return The assigned distributor.
    auto operator=(const CT3CDistributor& other) -> CT3CDistributor& = delete;

    /// @brief The default move assignment operator.
    /// @param other The distributor to be move assigned.
    /// @return The assigned distributor.
    auto operator=(CT3CDistributor&& other) noexcept -> CT3CDistributor& = delete;

    /// @brief The equality operator.
    /// @param other The distributor to be compared.
    /// @return True if distributors are equal, False otherwise.
    auto operator==(const CT3CDistributor& other) const -> bool = delete;

    /// @brief The equality operator.
    /// @param other The distributor to be compared.
    /// @return True if distributors are not equal, False otherwise.
    auto operator!=(const CT3CDistributor& other) const -> bool = delete;

    /// @brief Distributes buffer of integrals into local Fock matrix.
    /// @param buffer The integrals buffer.
    /// @param a_indices The compressed basis function indexes on center A.
    /// @param c_indices The compressed basis function indexes on center C.
    /// @param d_indices The compressed basis function indexes on center D.
    /// @param a_angmom The angular momentum of integrals buffer on center A.
    /// @param c_angmom The angular momentum of integrals buffer on center C.
    /// @param d_angmom Tthe angular momentum of integrals buffer on center D.
    /// @param ibra_gto The index of basis function on bra side.
    /// @param ket_range The index of the range [ket_first, ket_last) of basis functions on ket side.
    auto distribute(const CSimdArray<double>&        buffer,
                    const size_t                     offset,
                    const std::vector<size_t>&       a_indices,
                    const std::vector<size_t>&       c_indices,
                    const std::vector<size_t>&       d_indices,
                    const int                        a_angmom,
                    const int                        c_angmom,
                    const int                        d_angmom,
                    const size_t                     ibra_gto,
                    const std::pair<size_t, size_t>& ket_range) -> void;
    
   private:
    /// @brief The pointer to flat buffer.
    CT3FlatBuffer<double>* _t3_values;
};

#endif /* T3CDistributor_hpp */
