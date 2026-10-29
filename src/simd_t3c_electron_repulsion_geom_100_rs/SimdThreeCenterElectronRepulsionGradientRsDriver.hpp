//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//

#ifndef SimdThreeCenterElectronRepulsionGradientRsDriver_hpp
#define SimdThreeCenterElectronRepulsionGradientRsDriver_hpp

#include <cstddef>
#include <utility>
#include <vector>

#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "SparseTensor.hpp"
#include "TripleSparsityPattern.hpp"

/// @brief The derivative of the three-center electron repulsion integrals with
/// respect to the two atoms on bra side.
///
/// @note Six components an element: three Cartesian directions for each of the
/// two centers of the pair of basis functions. The center of the auxiliary
/// function is not differentiated, as the three derivatives of an integral sum
/// to zero and the caller takes the third from the other two.
///
/// @note The pattern is the one the calculation already holds. A gradient
/// screened on its own threshold would describe a different set of atom pairs
/// than the B vectors were formed over, and the two would no longer index the
/// same things. The note on the analytic gradient argues for a threshold two
/// decades tighter here; that is deliberately not done, and this is where it
/// would go.
class CSimdThreeCenterElectronRepulsionGradientRsDriver
{
   public:
    /// @brief Creates a driver with a target block size.
    /// @param block_size The target number of atom pairs of a block, or zero to
    /// choose it from the number of the threads and the number of the atom pairs.
    explicit CSimdThreeCenterElectronRepulsionGradientRsDriver(const size_t block_size = 0)
        : _block_size(block_size)
    {
    }

    CSimdThreeCenterElectronRepulsionGradientRsDriver(const CSimdThreeCenterElectronRepulsionGradientRsDriver &) = delete;
    CSimdThreeCenterElectronRepulsionGradientRsDriver(CSimdThreeCenterElectronRepulsionGradientRsDriver &&) noexcept = delete;
    ~CSimdThreeCenterElectronRepulsionGradientRsDriver() = default;
    auto operator=(const CSimdThreeCenterElectronRepulsionGradientRsDriver &)
        -> CSimdThreeCenterElectronRepulsionGradientRsDriver & = delete;
    auto operator=(CSimdThreeCenterElectronRepulsionGradientRsDriver &&) noexcept
        -> CSimdThreeCenterElectronRepulsionGradientRsDriver & = delete;

    /// @brief Gets target number of atom pairs of a block.
    auto
    get_block_size() const -> size_t
    {
        return _block_size;
    }

    /// @brief Computes the derivative of the integrals of the atoms on the
    /// auxiliary side, over a pattern the caller already holds.
    /// @param pattern The sparsity pattern, as the calculation formed it.
    /// @param molecule The molecule.
    /// @param basis The molecular basis on a and b sides.
    /// @param aux_basis The auxiliary molecular basis on c side.
    /// @return The sparse tensor of the derivative, of six components an element.
    /// @note The tensor is the size of the integrals times six. A caller which
    /// wants the whole auxiliary basis at once should think again: this is asked
    /// for one atom of the auxiliary side at a time and contracted at once.
    /// @param omega The range separation parameter, which must be positive.
    /// @return The derivatives of the Coulomb operator and of the attenuated one,
    /// as two sparse tensors on the one pattern, six components an element in each.
    /// @note They come back apart because the two carry different coefficients in a
    /// range separated functional and the caller is what knows them.
    /// @note Which three of the six blocks a kernel writes belong to which operator
    /// is not stated by the generator and cannot be read off the recursion. The
    /// constants below were **determined by comparison**: the Coulomb half must
    /// reproduce the unattenuated gradient driver, and only one assignment does.
    auto compute(const CTripleSparsityPattern &pattern,
                 const CMolecule              &molecule,
                 const CMolecularBasis        &basis,
                 const CMolecularBasis        &aux_basis,
                 const double                  omega) const -> std::pair<CSparseTensor, CSparseTensor>;

    /// @brief Gets the memory of the arena one thread holds while computing the
    /// derivative of the integrals of a sparsity pattern.
    /// @param pattern The sparsity pattern the derivative would be computed over.
    /// @param basis The molecular basis on a and b sides.
    /// @param aux_basis The auxiliary molecular basis on c side.
    /// @return The memory in bytes, which every thread of a call holds at once.
    /// @note The arena is the largest combination of the pattern by its widest
    /// block, as compute sizes it, so a caller adds this for every thread to the
    /// memory of the derivative itself to know the peak of a call.
    static auto arena_memory(const CTripleSparsityPattern &pattern,
                             const CMolecularBasis        &basis,
                             const CMolecularBasis        &aux_basis) -> size_t;

   private:
    /// @brief Gets the shape of the arena one thread holds for a sparsity pattern.
    /// @return The rows of the largest combination any block carries, of either
    /// differentiated center, and the atom pairs by the auxiliary atoms of the
    /// widest block.
    static auto _arena_shape(const CTripleSparsityPattern &pattern,
                             const CMolecularBasis        &basis,
                             const CMolecularBasis        &aux_basis) -> std::pair<size_t, size_t>;

    /// @brief Which of the six blocks a kernel writes belongs to which operator.
    static constexpr size_t _coulomb_block = 0;

    static constexpr size_t _attenuated_block = 3;

   private:
    /// @brief The target number of atom pairs of a block.
    size_t _block_size;
};

#endif /* SimdThreeCenterElectronRepulsionGradientRsDriver_hpp */
