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




#ifndef SimdT3CRsDistributor_hpp
#define SimdT3CRsDistributor_hpp

#include <algorithm>
#include <cstddef>
#include <vector>

#include <omp.h>

#include "AtomBasisTripleSparsity.hpp"
#include "OpenMPFunc.hpp"
#include "SparseTensor.hpp"

/// @brief Class CSimdT3CRsDistributor takes the integrals of a combination of basis
/// functions which a range separated three-center kernel writes as two blocks, and
/// puts each of them where it belongs.
///
/// The primary template does nothing; a consumer is served by a specialization.
///
/// @note Thread safety is a property of the specialization and is stated by each of
/// them.
template <class T>
class CSimdT3CRsDistributor
{
   public:
    /// @brief The constructor with the two storages to distribute into.
    CSimdT3CRsDistributor(T *coulomb, T *attenuated)

        : _coulomb{coulomb}

        , _attenuated{attenuated}
    {
    }

    CSimdT3CRsDistributor(const CSimdT3CRsDistributor &other) = delete;

    CSimdT3CRsDistributor(CSimdT3CRsDistributor &&other) noexcept = delete;

    ~CSimdT3CRsDistributor() = default;

    auto operator=(const CSimdT3CRsDistributor &other) -> CSimdT3CRsDistributor & = delete;

    auto operator=(CSimdT3CRsDistributor &&other) noexcept -> CSimdT3CRsDistributor & = delete;

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
    T *_coulomb;

    T *_attenuated;
};

/// @brief The specialization which distributes into a pair of sparse tensors.
///
/// @note Unlike the unattenuated distributor, this one cannot hand out a pointer
/// into its storage. A kernel writes the two operators as one contiguous block, the
/// Coulomb integrals first and the range separated ones after, while their
/// destinations are the same offset of two different tensors. So a scratch is handed
/// out and commit puts each half where it belongs. That copy is the price of holding
/// the two operators apart, and is the one thing this path pays which the
/// unattenuated one does not.
///
/// @note The scratch is one per thread, taken by the thread's own index, because
/// target and commit are called from inside the parallel loop over the combinations.
/// The outer vector is sized once at construction and never grown, so no thread
/// resizes what another is reading; each thread grows its own inner vector alone.
///
/// @note This specialization needs no synchronization of its own beyond that: the
/// blocks are disjoint in the values of each tensor, and so are the combinations of
/// basis functions within a block, exactly as in the unattenuated case.
template <>
class CSimdT3CRsDistributor<CSparseTensor>
{
   public:
    CSimdT3CRsDistributor(CSparseTensor *coulomb, CSparseTensor *attenuated)

        : _coulomb{coulomb}

        , _attenuated{attenuated}

        , _scratch(static_cast<size_t>(omp::get_number_of_threads()))
    {
    }

    CSimdT3CRsDistributor(const CSimdT3CRsDistributor &other) = delete;

    CSimdT3CRsDistributor(CSimdT3CRsDistributor &&other) noexcept = delete;

    ~CSimdT3CRsDistributor() = default;

    auto operator=(const CSimdT3CRsDistributor &other) -> CSimdT3CRsDistributor & = delete;

    auto operator=(CSimdT3CRsDistributor &&other) noexcept -> CSimdT3CRsDistributor & = delete;

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
        // NOTE: twice the values of the combination, as the kernel writes the two
        // operators one after the other into what it is given.
        const auto needed = 2 * ncomps * natoms * npairs;

        auto &mine = _scratch[thread()];

        if (mine.size() < needed) mine.resize(needed);

        return mine.data();
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
        const auto half = ncomps * natoms * npairs;

        const auto offset = block.element_offset(a_angular_momentum, a_index, b_angular_momentum, b_index,
                                                 c_angular_momentum, c_index);

        const auto *values = _scratch[thread()].data();

        std::copy_n(values, half, _coulomb->values(iblock) + offset);

        std::copy_n(values + half, half, _attenuated->values(iblock) + offset);
    }

   private:
    /// @brief The index of the calling thread, which is zero outside a parallel
    /// region and is where this thread's scratch lives.
    static auto
    thread() -> size_t
    {
        return static_cast<size_t>(omp_get_thread_num());
    }

    CSparseTensor *_coulomb;

    CSparseTensor *_attenuated;

    /// @brief One scratch per thread, each grown to the largest combination that
    /// thread has been handed.
    std::vector<std::vector<double>> _scratch;
};

#endif /* SimdT3CRsDistributor_hpp */
