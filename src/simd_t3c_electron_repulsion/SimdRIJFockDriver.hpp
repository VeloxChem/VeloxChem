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


#ifndef SimdRIJFockDriver_hpp
#define SimdRIJFockDriver_hpp

#include <cstddef>
#include <vector>

#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "PackedMatrix.hpp"
#include "SparseTensor.hpp"

/// @brief Class CSimdRIJFockDriver computes the B vectors of the resolution of the
/// identity approximation of the Coulomb interaction.
///
/// @note The B vectors are defined as B(q)_ij = sum over p of (ij|p) L(-1)_pq, where
/// (ij|p) are the three-center electron repulsion integrals and L(-1) is the inverse
/// of the two-center metric (p|q), which the caller supplies in the packed format.
///
/// @note The B vectors are screened on the atomic orbital pair, as the three-center
/// integrals are, and are dense in q, as the inverse metric is dense. Their memory is
/// therefore the surviving atomic orbital pairs times the auxiliary basis, which
/// reaches hundreds of gigabytes for a large molecule. The auxiliary functions of q
/// are what a calculation distributes over the MPI ranks, and the driver forms the
/// ones of the given atoms alone, so that a rank holds its own share and no more.
/// A single rank asks for all of them.
///
/// @note The sum over p runs over the whole auxiliary basis whatever share of q a
/// rank holds, so the three-center integrals cannot be distributed the same way. They
/// are formed in batches of the blocks of atomic orbital pairs instead, and each
/// batch is contracted and dropped before the next is formed, so that the peak is the
/// B vectors and one batch rather than the whole three-center tensor.
class CSimdRIJFockDriver
{
   public:
    /// @brief The default constructor.
    CSimdRIJFockDriver() = default;

    /// @brief Computes the B vectors of the resolution of the identity approximation.
    /// @param molecule The molecule to compute the B vectors of.
    /// @param basis The molecular basis on a and b sides.
    /// @param aux_basis The auxiliary molecular basis.
    /// @param inverse_metric The inverse of the two-center metric, in the packed
    /// format, as packlin::invert returns it.
    /// @param threshold The screening threshold.
    /// @param aux_atoms The atoms whose auxiliary functions the B vectors are formed
    /// for, as their indices in the molecule and without repetition, or empty for all
    /// of them.
    /// @return The sparse tensor of the B vectors, whose blocks hold the auxiliary
    /// functions of q in the place the three-center integrals hold those of p.
    auto compute_bq_vectors(const CMolecule        &molecule,
                            const CMolecularBasis  &basis,
                            const CMolecularBasis  &aux_basis,
                            const CPackedMatrix    &inverse_metric,
                            const double            threshold,
                            const std::vector<int> &aux_atoms = {}) const -> CSparseTensor;

   private:
    /// @brief The memory a batch of the three-center integrals is allowed to reach.
    /// @note The batch is the blocks of atomic orbital pairs whose integrals are
    /// formed at once. Making it larger costs memory and buys nothing beyond the
    /// point where the threads are busy, as the contraction of a batch is what is
    /// parallelized and not the batches.
    static constexpr size_t _batch_budget = size_t{4} * 1024 * 1024 * 1024;
};

#endif /* SimdRIJFockDriver_hpp */
