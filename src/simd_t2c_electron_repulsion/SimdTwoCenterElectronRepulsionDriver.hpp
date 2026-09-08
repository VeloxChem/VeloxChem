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



#ifndef SimdTwoCenterElectronRepulsionDriver_hpp
#define SimdTwoCenterElectronRepulsionDriver_hpp

#include <cstddef>
#include <vector>

#include "AtomBasisPairGroup.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "PackedMatrix.hpp"

/// @brief Class CSimdTwoCenterElectronRepulsionDriver computes the two-center
/// electron repulsion integrals of a molecular basis and stores them in a packed
/// matrix.
///
/// @note There is no sparsity pattern and no threshold. The Coulomb operator decays
/// as the inverse of the interatomic distance, so no atom pair of any molecule falls
/// below a threshold and the matrix is dense however large the molecule is. The
/// blocks of atom pairs are formed for the threads alone, by the geometry half of
/// sparsity::make_blocks, and none of them is described or screened.
///
/// @note There is no overload taking a pair of molecular bases, as the kinetic
/// energy driver has none. The repulsion of the charge distributions of two basis
/// functions of different molecular bases is not a quantity this driver is asked
/// for, unlike their overlap, so the driver takes one molecular basis alone.
class CSimdTwoCenterElectronRepulsionDriver
{
   public:
    /// @brief The constructor with target block size.
    /// @param block_size The target number of atom pairs of a block, or zero to
    /// choose it from the number of the threads and the number of the atom pairs.
    explicit CSimdTwoCenterElectronRepulsionDriver(const size_t block_size = 0)

        : _block_size(block_size)
    {
    }

    /// @brief Gets target number of atom pairs of a block.
    /// @return The target number of atom pairs, zero if it is chosen automatically.
    auto
    get_block_size() const -> size_t
    {
        return _block_size;
    }

    /// @brief Computes the two-center electron repulsion matrix of a molecular basis.
    /// @param molecule The molecule to compute the matrix of.
    /// @param basis The molecular basis on bra and ket sides.
    /// @return The symmetric packed matrix of the integrals.
    auto compute(const CMolecule &molecule, const CMolecularBasis &basis) const -> CPackedMatrix;

    /// @brief The largest target number of atom pairs of a block chosen when the
    /// size is not given.
    /// @note Fitted on fourteen threads over the six benchmark molecules in the def2
    /// universal jfit and jkfit sets, against the generated kernels. The block size
    /// is not what limits this driver: every ceiling from 64 to 512 lands within four
    /// to eight per cent of the best any of them reaches, where the wrong floor cost
    /// the overlap driver a factor of four. A ceiling of 128 wins the mean over the
    /// twelve cases and 256 wins their total time, because 128 buys 7 to 24 per cent
    /// on the three small molecules and costs 6 to 10 per cent on the three large
    /// ones, which carry the time. The value stays at 256 for that reason. Ceilings
    /// above 512 were not fitted: the curves are already rising there, the
    /// paracetamol cluster in jfit going from 35.7 ms at 512 to 41.7 at 1024 and
    /// ubiquitin from 481.7 to 526.4.
    /// @note The ceiling meets the floor. sparsity::min_block_size is 256 as well, so
    /// min(max(npairs / (blocks_per_thread * nthreads), 256), 256) is 256 for every
    /// molecule and the term which follows the size of the molecule never applies.
    /// The driver therefore divides tagrisso and ubiquitin into blocks of the same
    /// size, which the measurement says costs it little but which is a coincidence of
    /// the two constants rather than a choice. Refitting the floor moves this driver
    /// with it.
    static constexpr size_t max_block_size = 256;

   private:
    /// @brief Computes the integrals of the atom pairs of the blocks and adds them
    /// to the matrix.
    /// @param matrix The packed matrix to add the integrals to.
    /// @param molecule The molecule to compute the integrals of.
    /// @param basis The molecular basis on bra and ket sides.
    /// @param blocks The blocks of atom pairs to compute the integrals of.
    auto _compute_pair_blocks(CPackedMatrix                          &matrix,
                              const CMolecule                        &molecule,
                              const CMolecularBasis                  &basis,
                              const std::vector<CAtomBasisPairGroup> &blocks) const -> void;

    /// @brief Computes the integrals of the atoms with themselves and adds them to
    /// the matrix.
    /// @param matrix The packed matrix to add the integrals to.
    /// @param basis The molecular basis on bra and ket sides.
    /// @note The integral of two basis functions on the same atom does not depend on
    /// the position of the atom, so neither the molecule nor the blocks are needed
    /// here: the atoms of the molecular basis are enough.
    auto _compute_diagonal_blocks(CPackedMatrix &matrix, const CMolecularBasis &basis) const -> void;

    /// @brief The target number of atom pairs of a block, zero to choose it from
    /// the number of the threads and the number of the atom pairs.
    size_t _block_size;
};

#endif /* SimdTwoCenterElectronRepulsionDriver_hpp */
