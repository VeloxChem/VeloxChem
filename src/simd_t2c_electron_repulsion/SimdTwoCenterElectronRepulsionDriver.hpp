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
    /// @brief One combination of basis functions of one block, which is the unit of
    /// work the threads draw on.
    struct TPairTask
    {
        /// @brief The index of the block among the blocks of atom pairs.
        size_t iblock;

        /// @brief The index of the basis function on bra side within its atom basis.
        size_t i;

        /// @brief The index of the basis function on ket side within its atom basis.
        size_t j;
    };

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
    /// @note Fitted on fourteen threads over the benchmark molecules in the def2
    /// universal jfit and the correlation consistent rifit sets, against the kernels
    /// in the tree. A ceiling of 128 is the best of 64 to 1024 on the mean, on the
    /// worst case and on the total time, taking the worst case from 1.21 to 1.09. It
    /// helps the small and middle sized molecules, tagrisso in jfit going from 3.91
    /// ms to 3.24 and c60 in cc-pVQZ-rifit from 38.5 to 31.9, and costs the largest
    /// ones little. A ceiling of 64 measures the same, as the floor below never lets
    /// the size fall under 128 for any molecule here.
    /// @note The ceiling meets the floor. sparsity::min_block_size is 256, so
    /// min(max(npairs / (blocks_per_thread * nthreads), 256), 128) is 128 for every
    /// molecule and the term which follows the size of the molecule never applies.
    /// The driver divides tagrisso and ubiquitin into blocks of the same size, which
    /// the measurement says costs it little but which is a consequence of the two
    /// constants meeting rather than a choice. Refitting the floor moves this driver
    /// with it.
    /// @note What no choice of this constant reaches is that the best block size
    /// follows the basis as much as the molecule, and the formula sees only the count
    /// of atom pairs. On the screened path ubiquitin wants 16384 atom pairs a block
    /// in def2-svp and 2048 in cc-pV5Z, an eightfold difference on one molecule,
    /// because the cost of a pair differs by that much. Closing that would mean
    /// telling the block size something about the basis.
    static constexpr size_t max_block_size = 128;

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
