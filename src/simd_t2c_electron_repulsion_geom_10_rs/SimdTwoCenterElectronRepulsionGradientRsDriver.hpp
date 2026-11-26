//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//

#ifndef SimdTwoCenterElectronRepulsionGradientRsDriver_hpp
#define SimdTwoCenterElectronRepulsionGradientRsDriver_hpp

#include <cstddef>
#include <utility>
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
class CSimdTwoCenterElectronRepulsionGradientRsDriver
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
    explicit CSimdTwoCenterElectronRepulsionGradientRsDriver(const size_t block_size = 0)
        : _block_size(block_size)
    {
    }

    CSimdTwoCenterElectronRepulsionGradientRsDriver(const CSimdTwoCenterElectronRepulsionGradientRsDriver &) = delete;
    CSimdTwoCenterElectronRepulsionGradientRsDriver(CSimdTwoCenterElectronRepulsionGradientRsDriver &&) noexcept = delete;
    ~CSimdTwoCenterElectronRepulsionGradientRsDriver() = default;
    auto operator=(const CSimdTwoCenterElectronRepulsionGradientRsDriver &) -> CSimdTwoCenterElectronRepulsionGradientRsDriver & = delete;
    auto operator=(CSimdTwoCenterElectronRepulsionGradientRsDriver &&) noexcept -> CSimdTwoCenterElectronRepulsionGradientRsDriver & = delete;

    /// @brief Gets target number of atom pairs of a block.
    /// @return The target number of atom pairs, zero if it is chosen automatically.
    auto
    get_block_size() const -> size_t
    {
        return _block_size;
    }

    /// @brief Computes the gradients of the two-center integrals of the Coulomb
    /// operator and of the attenuated one, each contracted with its own weighted
    /// density.
    /// @param molecule The molecule.
    /// @param basis The molecular basis on bra and ket sides, which for the
    /// resolution of the identity is the auxiliary one.
    /// @param omega_coulomb The two-index fitted density of the Coulomb operator,
    /// symmetric and of the dimensions of the auxiliary basis.
    /// @param omega_attenuated The same for the attenuated operator. The two
    /// operators are fitted in metrics of their own and their weighted densities
    /// are different matrices; contracting either derivative with the other's
    /// would be a plausible looking gradient of nothing.
    /// @param omega The range separation parameter, which must be positive.
    /// @param atoms The atoms to compute the gradients of.
    /// @return The gradient of the Coulomb operator and that of the attenuated
    /// one, each a general matrix of one row of three components per atom, with the
    /// rows of the atoms not asked for left zero.
    /// @note They are returned apart rather than added, because the two carry
    /// different coefficients in a range separated functional and the caller is
    /// what knows them. Adding them here would need those coefficients here.
    /// @note What is returned is the sum over the pairs of auxiliary functions of
    /// the weighted density times the derivative of the integral, with no sign
    /// applied, as in the unattenuated driver.
    /// @note A pair of atoms is computed when **either** of its atoms is asked
    /// for, as the derivative of the pair contributes to both. The work therefore
    /// does not fall in proportion to the atoms asked for.
    auto compute(const CMolecule        &molecule,
                 const CMolecularBasis  &basis,
                 const CPackedMatrix    &omega_coulomb,
                 const CPackedMatrix    &omega_attenuated,
                 const double            omega,
                 const std::vector<int> &atoms) const -> std::pair<CPackedMatrix, CPackedMatrix>;

    /// @brief Computes them for every atom of the molecule.
    auto compute(const CMolecule       &molecule,
                 const CMolecularBasis &basis,
                 const CPackedMatrix   &omega_coulomb,
                 const CPackedMatrix   &omega_attenuated,
                 const double           omega) const -> std::pair<CPackedMatrix, CPackedMatrix>;

    /// @brief The largest target number of atom pairs of a block chosen when the
    /// block size is not named.
    static constexpr size_t max_block_size = 2048;

   private:
    /// @brief Adds the contribution of the blocks of atom pairs to the gradient.
    auto _compute_pair_blocks(CPackedMatrix                          &gradient_coulomb,
                              CPackedMatrix                          &gradient_attenuated,
                              const CMolecule                        &molecule,
                              const CMolecularBasis                  &basis,
                              const CPackedMatrix                    &omega_coulomb,
                              const CPackedMatrix                    &omega_attenuated,
                              const double                            omega,
                              const std::vector<CAtomBasisPairGroup> &blocks,
                              const std::vector<bool>                &wanted) const -> void;

    /// @brief Which of the six blocks a kernel writes belongs to which operator.
    /// @note The kernels write three Cartesian components of one operator and then
    /// three of the other, and which comes first is not stated by the generator and
    /// cannot be read off the recursion. These were **determined by comparison**:
    /// the Coulomb half must reproduce the unattenuated gradient driver exactly,
    /// and only one assignment does. Changing them on a hunch would be changing an
    /// answer that was measured.
    static constexpr size_t _coulomb_block = 0;

    static constexpr size_t _attenuated_block = 3;

    /// @brief The target number of atom pairs of a block.
    size_t _block_size;
};

#endif /* SimdTwoCenterElectronRepulsionGradientRsDriver_hpp */
