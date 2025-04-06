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

#ifndef T2CGeom10SumDistributor_hpp
#define T2CGeom10SumDistributor_hpp

#include <cstddef>
#include <utility>
#include <vector>

#include "SimdArray.hpp"

/// @brief Class CT3CGeom100SumDistributor provides methods for distributing vector of into flat buffer.
class CT2CGeom10SumDistributor
{
   public:
    /// @brief Creates an integral shells distributor.
    /// @param values  The pointer to gradient values.
    /// @param gamma The pointer to B_q vector data.
    CT2CGeom10SumDistributor(double* values, const double* gamma);

    /// @brief The default copy constructor.
    /// @param other The distributor to be copied.
    CT2CGeom10SumDistributor(const CT2CGeom10SumDistributor& other) = delete;

    /// @brief The default move constructor.
    /// @param other The distributor to be moved.
    CT2CGeom10SumDistributor(CT2CGeom10SumDistributor&& other) noexcept = delete;

    /// @brief The default destructor.
    ~CT2CGeom10SumDistributor() = default;

    /// @brief The default copy assignment operator.
    /// @param other The distributor to be copy assigned.
    /// @return The assigned distributor.
    auto operator=(const CT2CGeom10SumDistributor& other) -> CT2CGeom10SumDistributor& = delete;

    /// @brief The default move assignment operator.
    /// @param other The distributor to be move assigned.
    /// @return The assigned distributor.
    auto operator=(CT2CGeom10SumDistributor&& other) noexcept -> CT2CGeom10SumDistributor& = delete;

    /// @brief The equality operator.
    /// @param other The distributor to be compared.
    /// @return True if distributors are equal, False otherwise.
    auto operator==(const CT2CGeom10SumDistributor& other) const -> bool = delete;

    /// @brief The equality operator.
    /// @param other The distributor to be compared.
    /// @return True if distributors are not equal, False otherwise.
    auto operator!=(const CT2CGeom10SumDistributor& other) const -> bool = delete;
    
    /// @brief Distributes buffer of integrals into storage.
    /// @param buffer The integrals buffer.
    /// @param bra_indices The compressed contracted basis functions indexes on bra side.
    /// @param ket_indices The compressed contracted basis functions indexes on ket side.
    /// @param bra_angmom The angular momentum of integrals buffer on bra side.
    /// @param ket_angmom The angular momentum of integrals buffer on ket side.
    /// @param bra_igto The index of the basis function on bra side.
    /// @param ket_range The index of the range [ket_first, ket_last) of basis functions on ket side.
    /// @param diagonal True if basis functions blocks on bra and ket are the same, False otherwise.
    auto distribute(const CSimdArray<double>&        buffer,
                    const std::vector<size_t>&       bra_indices,
                    const std::vector<size_t>&       ket_indices,
                    const int                        bra_angmom,
                    const int                        ket_angmom,
                    const size_t                     bra_igto,
                    const std::pair<size_t, size_t>& ket_range,
                    const bool                       diagonal) -> void;
    
   private:
    /// @brief The pointer to gradient values.
    double* _grad_values;
    
    /// @brief The pointer to B_q vector values.
    const double* _ptr_gamma;
};

#endif /* T2CGeom10SumDistributor_hpp */
