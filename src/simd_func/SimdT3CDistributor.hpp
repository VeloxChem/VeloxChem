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



#ifndef SimdT3CDistributor_hpp
#define SimdT3CDistributor_hpp

#include <cstddef>

#include "AtomBasisTripleSparsity.hpp"
#include "SparseTensor.hpp"

/// @brief Class CSimdT3CDistributor hands a kernel the destination of the integrals
/// of one combination of basis functions of a three-center quantity and takes them
/// once they are written.
///
/// The destination is asked for rather than filled, so that a consumer whose layout
/// is the one the kernels write, such as a sparse tensor, is written into directly
/// and no copy is made. A consumer of a different shape, such as a contraction with
/// a fitting coefficient, hands out a scratch of its own and does its work in
/// commit, where the atom pairs and the atoms on c side of the block give the
/// mapping from the destination to the atoms.
///
/// The primary template does nothing; a consumer is served by a specialization.
///
/// @note Thread safety is a property of the specialization and is stated by each of
/// them. The blocks are handed out in parallel, so a specialization which
/// accumulates into storage shared between blocks carries its own reduction.
template <class T>
class CSimdT3CDistributor
{
   public:
    /// @brief The constructor with the storage to distribute into.
    /// @param storage The storage associated with the distributor.
    explicit CSimdT3CDistributor(T *storage)

        : _storage{storage}
    {
    }

    CSimdT3CDistributor(const CSimdT3CDistributor &other) = delete;

    CSimdT3CDistributor(CSimdT3CDistributor &&other) noexcept = delete;

    ~CSimdT3CDistributor() = default;

    auto operator=(const CSimdT3CDistributor &other) -> CSimdT3CDistributor & = delete;

    auto operator=(CSimdT3CDistributor &&other) noexcept -> CSimdT3CDistributor & = delete;

    /// @brief Gets destination of the integrals of a combination of basis functions
    /// of a block.
    /// @param block The sparsity pattern of the block.
    /// @param iblock The index of the block.
    /// @param a_angular_momentum The angular momentum of basis function on a side.
    /// @param a_index The index of basis function on a side.
    /// @param b_angular_momentum The angular momentum of basis function on b side.
    /// @param b_index The index of basis function on b side.
    /// @param c_angular_momentum The angular momentum of basis function on c side.
    /// @param c_index The index of basis function on c side.
    /// @param npairs The number of atom pairs surviving the screening of the
    /// combination.
    /// @param natoms The number of atoms on c side of the block.
    /// @param ncomps The number of triples of spherical components of the combination.
    /// @return The pointer to the values of the combination, with the atom pairs
    /// running fastest, the atom on c side next and the triple of angular components
    /// slowest.
    auto
    target(const CAtomBasisTripleSparsity &block,
           const size_t                    iblock,
           const int                       a_angular_momentum,
           const size_t                    a_index,
           const int                       b_angular_momentum,
           const size_t                    b_index,
           const int                       c_angular_momentum,
           const size_t                    c_index,
           const size_t                    npairs,
           const size_t                    natoms,
           const size_t                    ncomps) -> double *
    {
        return nullptr;
    }

    /// @brief Takes the integrals a kernel wrote into the destination above.
    /// @note The arguments are those of the matching call to target.
    auto
    commit(const CAtomBasisTripleSparsity &block,
           const size_t                    iblock,
           const int                       a_angular_momentum,
           const size_t                    a_index,
           const int                       b_angular_momentum,
           const size_t                    b_index,
           const int                       c_angular_momentum,
           const size_t                    c_index,
           const size_t                    npairs,
           const size_t                    natoms,
           const size_t                    ncomps) -> void
    {
    }

   protected:
    /// @brief The storage the distributor distributes into.
    T *_storage;
};

/// @brief The specialization which distributes into a sparse tensor.
/// @note The layout a kernel writes is the layout of the values of the tensor, so
/// the kernel writes into the tensor directly and commit has nothing to do.
/// @note This specialization needs no synchronization of its own: the blocks are
/// disjoint in the values of the tensor, and so are the combinations of basis
/// functions within a block.
template <>
class CSimdT3CDistributor<CSparseTensor>
{
   public:
    explicit CSimdT3CDistributor(CSparseTensor *storage)

        : _storage{storage}
    {
    }

    CSimdT3CDistributor(const CSimdT3CDistributor &other) = delete;

    CSimdT3CDistributor(CSimdT3CDistributor &&other) noexcept = delete;

    ~CSimdT3CDistributor() = default;

    auto operator=(const CSimdT3CDistributor &other) -> CSimdT3CDistributor & = delete;

    auto operator=(CSimdT3CDistributor &&other) noexcept -> CSimdT3CDistributor & = delete;

    auto
    target(const CAtomBasisTripleSparsity &block,
           const size_t                    iblock,
           const int                       a_angular_momentum,
           const size_t                    a_index,
           const int                       b_angular_momentum,
           const size_t                    b_index,
           const int                       c_angular_momentum,
           const size_t                    c_index,
           const size_t                    npairs,
           const size_t                    natoms,
           const size_t                    ncomps) -> double *
    {
        return _storage->values(iblock) + block.element_offset(a_angular_momentum, a_index, b_angular_momentum, b_index,
                                                               c_angular_momentum, c_index);
    }

    auto
    commit(const CAtomBasisTripleSparsity &block,
           const size_t                    iblock,
           const int                       a_angular_momentum,
           const size_t                    a_index,
           const int                       b_angular_momentum,
           const size_t                    b_index,
           const int                       c_angular_momentum,
           const size_t                    c_index,
           const size_t                    npairs,
           const size_t                    natoms,
           const size_t                    ncomps) -> void
    {
    }

   private:
    CSparseTensor *_storage;
};

#endif /* SimdT3CDistributor_hpp */
