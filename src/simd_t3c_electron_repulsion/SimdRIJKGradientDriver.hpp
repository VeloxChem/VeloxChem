//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//

#ifndef SimdRIJKGradientDriver_hpp
#define SimdRIJKGradientDriver_hpp

#include <cstddef>
#include <vector>

#include "DenseMatrix.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "PackedMatrix.hpp"
#include "SparseTensor.hpp"

/// @brief The Coulomb and exchange contributions to the molecular gradient,
/// through the resolution of the identity, from B vectors which are already
/// formed.
///
/// @note This driver computes nothing on its own that the Fock driver has
/// already computed. It is handed the B vectors and the metric that
/// CSimdRIJKFockDriver holds after a calculation, which is the expensive half of
/// the resolution of the identity and is the same for the energy and for its
/// gradient. What it forms are the derivative integrals and their contraction.
///
/// @note Threads and not ranks. The work divides over the atoms the gradient is
/// asked for and over the blocks of the sparsity pattern, both with OpenMP. The
/// arguments carry `aux_atoms` all the same: under MPI a rank holds the B vectors
/// of a share of the auxiliary atoms, and a gradient formed from that share is a
/// partial one which the ranks reduce, exactly as the Fock matrices are. Nothing
/// here divides anything over ranks; the argument is what lets a caller that does
/// say which share it holds.
class CSimdRIJKGradientDriver
{
   public:
    /// @brief Creates a gradient driver with a screening threshold and a target
    /// block size.
    /// @param threshold The screening threshold of the integrals.
    /// @param block_size The target number of atom pairs of a block, or zero to
    /// choose it from the number of the threads and the number of the atom pairs.
    explicit CSimdRIJKGradientDriver(const double threshold = 1.0e-12, const size_t block_size = 0)
        : _threshold(threshold)
        , _block_size(block_size)
    {
    }

    CSimdRIJKGradientDriver(const CSimdRIJKGradientDriver &other)                = delete;
    CSimdRIJKGradientDriver(CSimdRIJKGradientDriver &&other) noexcept            = delete;
    ~CSimdRIJKGradientDriver()                                                   = default;
    auto operator=(const CSimdRIJKGradientDriver &other) -> CSimdRIJKGradientDriver     & = delete;
    auto operator=(CSimdRIJKGradientDriver &&other) noexcept -> CSimdRIJKGradientDriver & = delete;
    auto operator==(const CSimdRIJKGradientDriver &other) const -> bool = delete;
    auto operator!=(const CSimdRIJKGradientDriver &other) const -> bool = delete;

    /// @brief Gets the screening threshold of the integrals.
    /// @return The screening threshold.
    auto get_threshold() const -> double;

    /// @brief Gets the target number of atom pairs of a block.
    /// @return The target number of atom pairs.
    auto get_block_size() const -> size_t;

    /// @brief Computes the Coulomb and exchange contributions to the gradient of
    /// the atoms it is asked for.
    /// @param molecule The molecule.
    /// @param basis The molecular basis.
    /// @param aux_basis The auxiliary molecular basis the identity is resolved in.
    /// @param bq_vectors The B vectors, as CSimdRIJKFockDriver formed them.
    /// @param metric The inverted Cholesky factor of the metric of the auxiliary
    /// basis, as that driver holds it.
    /// @param density The density matrix.
    /// @param coefficients The occupied molecular orbitals.
    /// @param atoms The atoms to compute the gradient of.
    /// @param aux_atoms The atoms of the auxiliary basis the B vectors span, or
    /// an empty list for all of them.
    /// @return The gradient, one row of three components per atom of the
    /// molecule, with the rows of the atoms not asked for left zero.
    /// @note The rows not asked for are zero rather than absent, so that a caller
    /// which asks for a share of the atoms can reduce what it is given without
    /// knowing which share the others took. A gradient of a subset does not sum
    /// to zero, and translational invariance is not a check on it.
    auto compute(const CMolecule        &molecule,
                 const CMolecularBasis  &basis,
                 const CMolecularBasis  &aux_basis,
                 const CSparseTensor    &bq_vectors,
                 const CPackedMatrix    &metric,
                 const CPackedMatrix    &density,
                 const CPackedMatrix    &coefficients,
                 const std::vector<int> &atoms,
                 const std::vector<int> &aux_atoms = {}) const -> CDenseMatrix;

    /// @brief Computes the Coulomb and exchange contributions to the gradient of
    /// every atom of the molecule.
    /// @return The gradient, one row of three components per atom.
    auto compute(const CMolecule       &molecule,
                 const CMolecularBasis &basis,
                 const CMolecularBasis &aux_basis,
                 const CSparseTensor   &bq_vectors,
                 const CPackedMatrix   &metric,
                 const CPackedMatrix   &density,
                 const CPackedMatrix   &coefficients) const -> CDenseMatrix;

   private:
    /// @brief The screening threshold of the integrals.
    double _threshold;

    /// @brief The target number of atom pairs of a block.
    size_t _block_size;
};

#endif /* SimdRIJKGradientDriver_hpp */
