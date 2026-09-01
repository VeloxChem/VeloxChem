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
/// @note The integrals are not screened and the driver carries no threshold. The
/// Coulomb operator decays as the inverse of the interatomic distance, so the
/// integral of an atom pair does not fall below a threshold at any distance a
/// molecule reaches: every atom pair and every combination of basis functions
/// contributes, the matrix is dense for every molecule, and the sparsity
/// patterns which the overlap and the kinetic energy are built on describe
/// nothing here. The atom basis pair groups are therefore divided into blocks
/// directly, and neither the ordering of the atom pairs by interatomic distance
/// nor the screening of the pairs of primitives is done.
/// @note The blocks divide the work alone, they do not divide the storage. Each
/// of them computes one combination of basis functions at a time into a buffer
/// of its own and adds it to the matrix, which the blocks share and write
/// disjoint elements of, as the atom pairs of two blocks are disjoint and the
/// basis functions of two atoms are as well.
class CSimdTwoCenterElectronRepulsionDriver
{
   public:
    /// @brief The number of blocks per thread aimed at when the target number of
    /// atom pairs of a block is chosen. It sets the size for the molecules which
    /// fall between the floor and the ceiling below, and those alone: a small
    /// molecule is held at the floor and a large one at the ceiling whatever this
    /// is, so it decides little here.
    static constexpr size_t blocks_per_thread = 2;

    /// @brief The smallest target number of atom pairs of a block chosen. A block
    /// carries a fixed cost which does not shrink with the atom pairs it holds,
    /// chiefly the coordinates and the solid harmonics it forms, so a molecule
    /// too small to fill the threads is divided into fewer blocks rather than
    /// into blocks whose fixed cost outweighs their work.
    /// @note This is far below the two thousand and forty eight of the sparse
    /// matrix, whose blocks carry the bisection of the screening as their fixed
    /// cost. Nothing is screened here and a block does much more arithmetic per
    /// atom pair, so its fixed cost is repaid by far fewer of them. Measured on
    /// fourteen threads over the def2 universal fitting sets, it is worth twice
    /// the throughput of two thousand and forty eight on tagrisso and taxol.
    /// Sixty four is the better choice below about forty atoms, by a third on a
    /// molecule of thirty, and the worse one from there up, by three to ten per
    /// cent; the molecules it wins are already under two milliseconds. It binds
    /// for the small molecules alone, as a large one asks for blocks larger than
    /// this of its own accord.
    static constexpr size_t min_block_size = 128;

    /// @brief The largest target number of atom pairs of a block chosen. The size
    /// computed from the threads grows with the molecule, which suits the sparse
    /// matrix, whose blocks carry the bisection of the screening and are worth
    /// making large. A block here is cheap to start, so a large molecule is
    /// better divided into many small blocks than into a few per thread: the
    /// dynamic loop then has enough of them to even out the ones which differ in
    /// cost, and each of them keeps its Boys function and its buffer in the cache.
    /// @note Without this the size followed the threads alone and reached 7348
    /// atom pairs for crambin and 27038 for ubiquitin, one to two orders of
    /// magnitude above what the timings want. Measured on fourteen threads over
    /// the def2 universal fitting sets, the time is flat below about three
    /// hundred and rises steeply above it: crambin in jkfit takes 271 ms at 7348
    /// atom pairs a block, 195 at 918, 159 at 306 and 158 at 153, and ubiquitin
    /// in jkfit 1253 ms at 27038 and 731 at 211. The value is taken inside the
    /// flat range rather than at the floor, so that a molecule between the two is
    /// not divided more finely than it repays.
    static constexpr size_t max_block_size = 256;

    /// @brief Computes the two-center electron repulsion matrix of a molecular
    /// basis.
    /// @param molecule The molecule to compute the matrix of.
    /// @param basis The molecular basis on bra and ket sides.
    /// @return The symmetric packed matrix of the integrals.
    auto compute(const CMolecule &molecule, const CMolecularBasis &basis) const -> CPackedMatrix;

    // NOTE: there is no overload taking a pair of molecular bases, as the
    // kinetic energy driver has none. The repulsion of the charge distributions
    // of two basis functions of different molecular bases is not a quantity this
    // driver is asked for, unlike their overlap, so the driver takes one
    // molecular basis alone.

   private:
    /// @brief Computes the integrals of the atom pairs of the blocks and adds
    /// them to the matrix.
    /// @param matrix The packed matrix to add the integrals to.
    /// @param molecule The molecule to compute the integrals of.
    /// @param basis The molecular basis on bra and ket sides.
    /// @param blocks The blocks of atom pairs to compute the integrals of.
    auto _compute_pair_blocks(CPackedMatrix                          &matrix,
                              const CMolecule                        &molecule,
                              const CMolecularBasis                  &basis,
                              const std::vector<CAtomBasisPairGroup> &blocks) const -> void;

    /// @brief Computes the integrals of the atoms of the diagonal atom pairs of
    /// the blocks and adds them to the matrix.
    /// @param matrix The packed matrix to add the integrals to.
    /// @param basis The molecular basis on bra and ket sides.
    /// @param blocks The blocks of atom pairs to compute the integrals of.
    auto _compute_diagonal_blocks(CPackedMatrix                          &matrix,
                                  const CMolecularBasis                  &basis,
                                  const std::vector<CAtomBasisPairGroup> &blocks) const -> void;
};

#endif /* SimdTwoCenterElectronRepulsionDriver_hpp */
