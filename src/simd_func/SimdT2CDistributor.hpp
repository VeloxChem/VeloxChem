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


#ifndef SimdT2CDistributor_hpp
#define SimdT2CDistributor_hpp

#include <cstddef>

#include "AtomBasisDiagonalSparsity.hpp"
#include "AtomBasisPairSparsity.hpp"
#include "SparseMatrix.hpp"

/// @brief Class CSimdT2CDistributor hands a kernel the destination of the
/// integrals of one combination of basis functions and takes them once they are
/// written.
///
/// The destination is asked for rather than filled, so that a consumer whose
/// layout is the one the kernels write, such as a sparse matrix, is written into
/// directly and no copy is made. A consumer of a different shape, such as a
/// contraction with a density or with response vectors, hands out a scratch of its
/// own and does its work in commit, where the atom pairs of the block give the
/// mapping from the columns of the destination to the atoms.
///
/// The diagonal blocks are served by the same pair of calls, so that such a
/// consumer is notified of the integrals of the atom pairs of an atom with itself
/// as it is of the integrals of the other atom pairs.
///
/// The primary template does nothing; a consumer is served by a specialization.
///
/// @note Thread safety is a property of the specialization and is stated by each
/// of them. The blocks are handed out in parallel, so a specialization which
/// accumulates into storage shared between blocks carries its own reduction.
template <class T>
class CSimdT2CDistributor
{
   public:
    /// @brief The constructor with the storage to distribute into.
    /// @param storage The storage associated with the distributor.
    explicit CSimdT2CDistributor(T *storage)

        : _storage{storage}
    {
    }

    CSimdT2CDistributor(const CSimdT2CDistributor &other) = delete;

    CSimdT2CDistributor(CSimdT2CDistributor &&other) noexcept = delete;

    ~CSimdT2CDistributor() = default;

    auto operator=(const CSimdT2CDistributor &other) -> CSimdT2CDistributor & = delete;

    auto operator=(CSimdT2CDistributor &&other) noexcept -> CSimdT2CDistributor & = delete;

    /// @brief Gets destination of the integrals of a combination of basis functions
    /// of an off-diagonal block.
    /// @param block The sparsity pattern of the off-diagonal block.
    /// @param iblock The index of the off-diagonal block.
    /// @param bra_angular_momentum The angular momentum of basis function on bra side.
    /// @param bra_index The index of basis function on bra side.
    /// @param ket_angular_momentum The angular momentum of basis function on ket side.
    /// @param ket_index The index of basis function on ket side.
    /// @param nvalues The number of atom pairs surviving the screening of the combination.
    /// @param ncomps The number of spherical components of the combination.
    /// @return The pointer to ncomps rows of nvalues columns, the components of the
    /// bra side running slowest.
    auto
    target(const CAtomBasisPairSparsity &block,
           const size_t                  iblock,
           const int                     bra_angular_momentum,
           const size_t                  bra_index,
           const int                     ket_angular_momentum,
           const size_t                  ket_index,
           const size_t                  nvalues,
           const size_t                  ncomps) -> double *
    {
        return nullptr;
    }

    /// @brief Takes the integrals a kernel wrote into the destination above.
    /// @note The arguments are those of the matching call to target.
    auto
    commit(const CAtomBasisPairSparsity &block,
           const size_t                  iblock,
           const int                     bra_angular_momentum,
           const size_t                  bra_index,
           const int                     ket_angular_momentum,
           const size_t                  ket_index,
           const size_t                  nvalues,
           const size_t                  ncomps) -> void
    {
    }

    /// @brief Gets destination of the integral of a combination of basis functions
    /// of a diagonal block, which is a single value.
    /// @param block The sparsity pattern of the diagonal block.
    /// @param iblock The index of the diagonal block.
    /// @param bra_angular_momentum The angular momentum of basis function on bra side.
    /// @param bra_index The index of basis function on bra side.
    /// @param ket_angular_momentum The angular momentum of basis function on ket side.
    /// @param ket_index The index of basis function on ket side.
    /// @return The pointer to the value.
    /// @note The overlap of two basis functions on the same atom does not depend on
    /// the position of the atom, so the atoms of the block carry the value rather
    /// than a value being stored for each of them.
    auto
    diagonal_target(const CAtomBasisDiagonalSparsity &block,
                    const size_t                      iblock,
                    const int                         bra_angular_momentum,
                    const size_t                      bra_index,
                    const int                         ket_angular_momentum,
                    const size_t                      ket_index) -> double *
    {
        return nullptr;
    }

    /// @brief Takes the integral written into the destination above.
    /// @note The arguments are those of the matching call to diagonal_target.
    auto
    diagonal_commit(const CAtomBasisDiagonalSparsity &block,
                    const size_t                      iblock,
                    const int                         bra_angular_momentum,
                    const size_t                      bra_index,
                    const int                         ket_angular_momentum,
                    const size_t                      ket_index) -> void
    {
    }

   protected:
    /// @brief The storage the distributor distributes into.
    T *_storage;
};

/// @brief The specialization which distributes into a sparse matrix.
/// @note The layout a kernel writes is the layout of the values of the matrix, so
/// the kernel writes into the matrix directly and neither commit has anything to do.
/// @note This specialization needs no synchronization of its own: the blocks are
/// disjoint in the values of the matrix, and so are the combinations of basis
/// functions within a block.
template <>
class CSimdT2CDistributor<CSparseMatrix>
{
   public:
    explicit CSimdT2CDistributor(CSparseMatrix *storage)

        : _storage{storage}
    {
    }

    CSimdT2CDistributor(const CSimdT2CDistributor &other) = delete;

    CSimdT2CDistributor(CSimdT2CDistributor &&other) noexcept = delete;

    ~CSimdT2CDistributor() = default;

    auto operator=(const CSimdT2CDistributor &other) -> CSimdT2CDistributor & = delete;

    auto operator=(CSimdT2CDistributor &&other) noexcept -> CSimdT2CDistributor & = delete;

    auto
    target(const CAtomBasisPairSparsity &block,
           const size_t                  iblock,
           const int                     bra_angular_momentum,
           const size_t                  bra_index,
           const int                     ket_angular_momentum,
           const size_t                  ket_index,
           const size_t                  nvalues,
           const size_t                  ncomps) -> double *
    {
        return _storage->pair_values(iblock, bra_angular_momentum, bra_index, ket_angular_momentum, ket_index);
    }

    auto
    commit(const CAtomBasisPairSparsity &block,
           const size_t                  iblock,
           const int                     bra_angular_momentum,
           const size_t                  bra_index,
           const int                     ket_angular_momentum,
           const size_t                  ket_index,
           const size_t                  nvalues,
           const size_t                  ncomps) -> void
    {
    }

    auto
    diagonal_target(const CAtomBasisDiagonalSparsity &block,
                    const size_t                      iblock,
                    const int                         bra_angular_momentum,
                    const size_t                      bra_index,
                    const int                         ket_angular_momentum,
                    const size_t                      ket_index) -> double *
    {
        return _storage->diagonal_values(iblock, bra_angular_momentum, bra_index, ket_angular_momentum, ket_index);
    }

    auto
    diagonal_commit(const CAtomBasisDiagonalSparsity &block,
                    const size_t                      iblock,
                    const int                         bra_angular_momentum,
                    const size_t                      bra_index,
                    const int                         ket_angular_momentum,
                    const size_t                      ket_index) -> void
    {
    }

   private:
    CSparseMatrix *_storage;
};

#endif /* SimdT2CDistributor_hpp */
