//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//

#ifndef SimdTwoCenterElectronRepulsionGradientDriver_hpp
#define SimdTwoCenterElectronRepulsionGradientDriver_hpp

#include <cstddef>
#include <vector>

#include "AtomBasisPairGroup.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "PackedMatrix.hpp"

/// @brief The gradient of the two-center electron repulsion integrals of an
/// auxiliary basis, contracted against Omega as it is formed.
///
/// @note The derivative is never stored. A gradient of the two-center integrals
/// would be three numbers for every pair of auxiliary functions, which is three
/// times the metric; the caller wants three numbers per atom, so each batch of
/// derivatives is reduced against Omega and discarded.
///
/// @note Nothing is screened, which is the same choice the driver of the
/// integrals themselves makes: the Coulomb operator falls off as one over the
/// distance and no pair of atoms of it is negligible. The derivative falls off
/// one power faster and is no sparser in any way that pays for a pattern.
///
/// @note The derivative of the second center is not computed. The two centers of
/// a two-center integral move against each other, so their derivatives sum to
/// zero, and a pair is swept once: what is formed for the atom on bra side is
/// subtracted from the atom on ket side.
class CSimdTwoCenterElectronRepulsionGradientDriver
{
    /// @brief One combination of basis functions of one block, which is the unit
    /// of work the threads draw on.
    struct TPairTask
    {
        size_t iblock;

        size_t i;

        size_t j;
    };

   public:
    /// @brief The constructor with target block size.
    /// @param block_size The target number of atom pairs of a block, or zero to
    /// choose it from the number of the threads and the number of the atom pairs.
    explicit CSimdTwoCenterElectronRepulsionGradientDriver(const size_t block_size = 0)
        : _block_size(block_size)
    {
    }

    CSimdTwoCenterElectronRepulsionGradientDriver(const CSimdTwoCenterElectronRepulsionGradientDriver &) = delete;
    CSimdTwoCenterElectronRepulsionGradientDriver(CSimdTwoCenterElectronRepulsionGradientDriver &&) noexcept = delete;
    ~CSimdTwoCenterElectronRepulsionGradientDriver() = default;
    auto operator=(const CSimdTwoCenterElectronRepulsionGradientDriver &) -> CSimdTwoCenterElectronRepulsionGradientDriver & = delete;
    auto operator=(CSimdTwoCenterElectronRepulsionGradientDriver &&) noexcept -> CSimdTwoCenterElectronRepulsionGradientDriver & = delete;

    /// @brief Gets target number of atom pairs of a block.
    /// @return The target number of atom pairs, zero if it is chosen automatically.
    auto
    get_block_size() const -> size_t
    {
        return _block_size;
    }

    /// @brief Computes the gradient of the two-center electron repulsion integrals
    /// contracted with Omega.
    /// @param molecule The molecule.
    /// @param basis The molecular basis on bra and ket sides, which for the
    /// resolution of the identity is the auxiliary one.
    /// @param omega The two-index fitted density, symmetric and of the dimensions
    /// of the auxiliary basis, as the first phase of the gradient forms it.
    /// @param atoms The atoms to compute the gradient of.
    /// @return The gradient, a general matrix of one row of three components per
    /// atom of the molecule, with the rows of the atoms not asked for left zero.
    /// @note What is returned is the sum over the pairs of auxiliary functions of
    /// Omega times the derivative of the integral, with no sign applied. The term
    /// of an RI-JK gradient is the negative of it, and the caller which knows
    /// which term it is asking for applies the sign.
    /// @note A pair of atoms is computed when **either** of its atoms is asked
    /// for, as the derivative of the pair contributes to both. The work therefore
    /// does not fall in proportion to the atoms asked for.
    auto compute(const CMolecule        &molecule,
                 const CMolecularBasis  &basis,
                 const CPackedMatrix    &omega,
                 const std::vector<int> &atoms) const -> CPackedMatrix;

    /// @brief Computes it for every atom of the molecule.
    auto compute(const CMolecule &molecule, const CMolecularBasis &basis, const CPackedMatrix &omega) const
        -> CPackedMatrix;

    /// @brief The largest target number of atom pairs of a block chosen when the
    /// block size is not named.
    static constexpr size_t max_block_size = 2048;

   private:
    /// @brief Adds the contribution of the blocks of atom pairs to the gradient.
    auto _compute_pair_blocks(CPackedMatrix                          &gradient,
                              const CMolecule                        &molecule,
                              const CMolecularBasis                  &basis,
                              const CPackedMatrix                    &omega,
                              const std::vector<CAtomBasisPairGroup> &blocks,
                              const std::vector<bool>                &wanted) const -> void;

    /// @brief The target number of atom pairs of a block.
    size_t _block_size;
};

#endif /* SimdTwoCenterElectronRepulsionGradientDriver_hpp */
