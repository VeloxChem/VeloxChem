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

#ifndef T4CEriTensorDistributor_hpp
#define T4CEriTensorDistributor_hpp

#include <array>
#include <string>
#include <vector>

#include "Dense4DTensor.hpp"
#include "GtoPairBlock.hpp"
#include "Matrices.hpp"
#include "Matrix.hpp"
#include "SimdArray.hpp"

/// @brief Class CT4CEriTensorDistributor provides methods for distributing single Fock matrix associated
/// with single density matrix.
class CT4CEriTensorDistributor
{
   public:
    /// @brief The default destructor.
    CT4CEriTensorDistributor() = default;

    /// Creates a Fock matrix distributor.
    /// @param eri_tensor The ERI tensor.
    CT4CEriTensorDistributor(CDense4DTensor* eri_tensor);

    /// @brief The default copy constructor.
    /// @param other The distributor to be copied.
    CT4CEriTensorDistributor(const CT4CEriTensorDistributor& other) = delete;

    /// @brief The default move constructor.
    /// @param other The distributor to be moved.
    CT4CEriTensorDistributor(CT4CEriTensorDistributor&& other) noexcept = delete;

    /// @brief The default destructor.
    ~CT4CEriTensorDistributor() = default;

    /// @brief The default copy assignment operator.
    /// @param other The distributor to be copy assigned.
    /// @return The assigned distributor.
    auto operator=(const CT4CEriTensorDistributor& other) -> CT4CEriTensorDistributor& = delete;

    /// @brief The default move assignment operator.
    /// @param other The distributor to be move assigned.
    /// @return The assigned distributor.
    auto operator=(CT4CEriTensorDistributor&& other) noexcept -> CT4CEriTensorDistributor& = delete;

    /// @brief The equality operator.
    /// @param other The distributor to be compared.
    /// @return True if distributors are equal, False otherwise.
    auto operator==(const CT4CEriTensorDistributor& other) const -> bool = delete;

    /// @brief The equality operator.
    /// @param other The distributor to be compared.
    /// @return True if distributors are not equal, False otherwise.
    auto operator!=(const CT4CEriTensorDistributor& other) const -> bool = delete;

    /// Sets local matrices and their local/global indices.
    /// @param bra_gto_pair_block The basis function pairs block on bra side.
    /// @param ket_gto_pair_block The basis function pairs block on ket side.
    auto set_indices(const CGtoPairBlock& bra_gto_pair_block, const CGtoPairBlock& ket_gto_pair_block) -> void;

    /// @brief Distributes buffer of integrals into local Fock matrix.
    /// @param buffer The integrals buffer.
    /// @param a_indices The compressed basis function indexes on center A.
    /// @param b_indices The compressed basis function indexes on center B.
    /// @param c_indices The compressed basis function indexes on center C.
    /// @param d_indices The compressed basis function indexes on center D.
    /// @param a_angmom The angular momentum of integrals buffer on center A.
    /// @param b_angmom The angular momentum of integrals buffer on center B.
    /// @param c_angmom The angular momentum of integrals buffer on center C.
    /// @param d_angmom Tthe angular momentum of integrals buffer on center D.
    /// @param ibra_gto The index of basis function on bra side.
    /// @param ket_range The index of the range [ket_first, ket_last) of basis functions on ket side.
    /// @param diagonal The flag to indicate diagonal subblock of matrix.
    auto distribute(const CSimdArray<double>&        buffer,
                    const size_t                     offset,
                    const std::vector<size_t>&       a_indices,
                    const std::vector<size_t>&       b_indices,
                    const std::vector<size_t>&       c_indices,
                    const std::vector<size_t>&       d_indices,
                    const int                        a_angmom,
                    const int                        b_angmom,
                    const int                        c_angmom,
                    const int                        d_angmom,
                    const size_t                     ibra_gto,
                    const std::pair<size_t, size_t>& ket_range,
                    const bool                       diagonal) -> void;

    auto need_omega() const -> bool;
    auto get_omega() const -> double;

   private:
    /// @brief The Fock matrix associated with distributor.
    CDense4DTensor* _eri_tensor;

    /// @brief The local indices for center A.
    std::vector<size_t> _a_loc_indices;

    /// @brief The local indices for center B.
    std::vector<size_t> _b_loc_indices;

    /// @brief The local indices for center C.
    std::vector<size_t> _c_loc_indices;

    /// @brief The local indices for center D.
    std::vector<size_t> _d_loc_indices;

    /// @brief The global indices for center A.
    std::vector<size_t> _a_glob_indices;

    /// @brief The global indices for center B.
    std::vector<size_t> _b_glob_indices;

    /// @brief The global indices for center C.
    std::vector<size_t> _c_glob_indices;

    /// @brief The global indices for center D.
    std::vector<size_t> _d_glob_indices;
};

#endif /* T4CEriTensorDistributor_hpp */
