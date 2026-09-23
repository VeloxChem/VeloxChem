//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//

#ifndef SimdRIJGradientDriver_hpp
#define SimdRIJGradientDriver_hpp

#include <cstddef>
#include <vector>

#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "PackedMatrix.hpp"

/// @brief The Coulomb contribution to the molecular gradient, through the
/// resolution of the identity, for a functional which asks for no exact exchange.
///
/// @note What this driver needs of a calculation is one vector: the fitting
/// coefficients gamma, which CSimdRIJFockDriver already solved for on every
/// iteration and which the converged density fixes. There are no B vectors to be
/// handed, no metric, and no occupied orbitals -- the fitting of a Coulomb matrix
/// is closed over the basis functions and never over an orbital, so nothing of the
/// wavefunction beyond the density enters. The driver which also fits the exchange
/// has to be handed all four.
///
/// @note The gradient is two terms:
///
///     dE/dx = 4 sum over P of gamma(P) sum over mn of D(mn) d(mn|P)/dx
///           - 2 sum over PQ of gamma(P) gamma(Q) dJ(PQ)/dx
///
/// the first a derivative of the three-center integrals contracted against the
/// density and the coefficients, the second the derivative of the metric against
/// their outer product. Those factors are the closed shell's, where the density
/// handed in is one spin's and the coefficients are of that density.
///
/// @note An open shell fits the **total** density and carries one and a half in
/// place of the four and the two, which `compute_open_shell` applies. Setting the
/// two spins equal returns the closed shell expression exactly: the total density
/// is twice one spin's and so are its coefficients, and 4 times a half times a
/// half is one, as 2 times a quarter is a half.
///
/// @note Threads and not ranks. The work divides over the atoms of the auxiliary
/// basis and over the blocks of the sparsity pattern. `aux_atoms` carries the share
/// a rank holds: a gradient formed from a share is a partial one which the ranks
/// reduce, exactly as the Coulomb matrices are, and nothing here divides anything
/// over ranks by itself.
class CSimdRIJGradientDriver
{
   public:
    /// @brief Creates a gradient driver with a screening threshold and a target
    /// block size.
    /// @param threshold The screening threshold of the integrals, which must be
    /// the one the calculation formed its fitting with, so that the derivative
    /// runs over the atom pairs the energy did.
    /// @param block_size The target number of atom pairs of a block, or zero to
    /// choose it from the number of the threads and the number of the atom pairs.
    CSimdRIJGradientDriver(const double threshold = 1.0e-12, const size_t block_size = 0)
        : _threshold(threshold)
        , _block_size(block_size)
    {
    }

    CSimdRIJGradientDriver(const CSimdRIJGradientDriver &other)                          = delete;
    CSimdRIJGradientDriver(CSimdRIJGradientDriver &&other) noexcept                      = delete;
    ~CSimdRIJGradientDriver()                                                            = default;
    auto operator=(const CSimdRIJGradientDriver &other) -> CSimdRIJGradientDriver      & = delete;
    auto operator=(CSimdRIJGradientDriver &&other) noexcept -> CSimdRIJGradientDriver  & = delete;
    auto operator==(const CSimdRIJGradientDriver &other) const -> bool                   = delete;
    auto operator!=(const CSimdRIJGradientDriver &other) const -> bool                   = delete;

    /// @brief Gets the screening threshold of the integrals.
    /// @return The screening threshold.
    auto get_threshold() const -> double;

    /// @brief Gets the target number of atom pairs of a block.
    /// @return The target number of atom pairs.
    auto get_block_size() const -> size_t;

    /// @brief Computes the Coulomb contribution to the gradient of the given atoms.
    /// @param molecule The molecule.
    /// @param basis The molecular basis.
    /// @param aux_basis The auxiliary molecular basis the identity is resolved in.
    /// @param fitting The fitting coefficients, gamma, as CSimdRIJFockDriver's
    /// solve_fitting answers them for the converged density.
    /// @param density The density matrix of one spin, which is the density the
    /// fitting coefficients were solved from.
    /// @param atoms The atoms to compute the gradient of.
    /// @param aux_atoms The atoms of the auxiliary basis this rank holds. An empty
    /// list is an empty share and not every atom: a caller which wants every atom
    /// says so, or takes the overload which does it for them. See the note in
    /// SimdRIJGradientDriver::_three_center for what reading it the other way cost.
    /// @param with_metric Whether to add the two-center term.
    /// @return The gradient, a general matrix of one row of three components per
    /// atom of the molecule, with the rows of the atoms not asked for left zero.
    /// @note The rows not asked for are zero rather than absent, so that a caller
    /// which asks for a share of the atoms can add what it gets to a whole
    /// gradient without knowing where its share sat.
    /// @note `with_metric` exists because the fitting coefficients are **not**
    /// divided between the ranks the way the B vectors of an exchange are: the
    /// fitting is solved after the ranks have added their halves of its right hand
    /// side, so every rank holds all of it and would compute the whole two-center
    /// term. Ranks which reduce their gradients afterwards must therefore have
    /// exactly one of them ask for it. The three-center term needs no such flag:
    /// `aux_atoms` already says which share of it a rank is forming.
    auto compute(const CMolecule           &molecule,
                 const CMolecularBasis     &basis,
                 const CMolecularBasis     &aux_basis,
                 const std::vector<double> &fitting,
                 const CPackedMatrix       &density,
                 const std::vector<int>    &atoms,
                 const std::vector<int>    &aux_atoms  = {},
                 const bool                 with_metric = true) const -> CPackedMatrix;

    /// @brief Computes the Coulomb contribution to the gradient of every atom.
    /// @return The gradient, a general matrix of three components per atom.
    auto compute(const CMolecule           &molecule,
                 const CMolecularBasis     &basis,
                 const CMolecularBasis     &aux_basis,
                 const std::vector<double> &fitting,
                 const CPackedMatrix       &density) const -> CPackedMatrix;

    /// @brief Computes the Coulomb contribution to the gradient of an open shell.
    /// @param density The **total** density, of both spins added, which is the
    /// density an open shell fits and the one its coefficients are of.
    /// @param fitting The fitting coefficients of that total density.
    /// @param atoms The atoms to compute the gradient of.
    /// @param aux_atoms The atoms of the auxiliary basis this rank holds, an empty
    /// list being an empty share.
    /// @param with_metric Whether to add the two-center term.
    /// @return The gradient, one row per atom.
    /// @note A separate routine and not a flag on the one above. The two differ in
    /// which density they are handed as much as in the factors they carry, and a
    /// caller which passed one spin's density here would be told nothing: the
    /// answer would simply be a quarter of the right one.
    auto compute_open_shell(const CMolecule           &molecule,
                            const CMolecularBasis     &basis,
                            const CMolecularBasis     &aux_basis,
                            const std::vector<double> &fitting,
                            const CPackedMatrix       &density,
                            const std::vector<int>    &atoms,
                            const std::vector<int>    &aux_atoms  = {},
                            const bool                 with_metric = true) const -> CPackedMatrix;

    /// @brief Computes the Coulomb contribution to an open shell's gradient over
    /// every atom of the molecule.
    auto compute_open_shell(const CMolecule           &molecule,
                            const CMolecularBasis     &basis,
                            const CMolecularBasis     &aux_basis,
                            const std::vector<double> &fitting,
                            const CPackedMatrix       &density) const -> CPackedMatrix;

   private:
    /// @brief Adds the three-center term to the gradient, one atom of the
    /// auxiliary basis at a time.
    /// @note The derivative integrals of an auxiliary atom are formed, contracted
    /// and discarded before the next atom is asked for. Holding the whole
    /// auxiliary basis at once is the basis functions squared times the auxiliary
    /// basis times six, which is hundreds of gigabytes on anything of a size worth
    /// fitting.
    /// @param coulomb_factor Four for a closed shell and one for an open one.
    auto _three_center(CPackedMatrix             &gradient,
                       const CMolecule           &molecule,
                       const CMolecularBasis     &basis,
                       const CMolecularBasis     &aux_basis,
                       const std::vector<double> &fitting,
                       const CPackedMatrix       &density,
                       const double               coulomb_factor,
                       const std::vector<bool>   &wanted,
                       const std::vector<int>    &aux_atoms) const -> void;

    /// @brief The two-index fitted density, the outer product of the fitting
    /// coefficients which the derivative of the metric is contracted against.
    /// @param fitting The fitting coefficients.
    /// @param factor Two for a closed shell and a half for an open one.
    /// @return Omega, symmetric and of the dimensions of the auxiliary basis.
    auto _omega(const std::vector<double> &fitting, const double factor) const -> CPackedMatrix;

    /// @brief The body both spin cases share, which differ only in their factors.
    auto _compute(const CMolecule           &molecule,
                  const CMolecularBasis     &basis,
                  const CMolecularBasis     &aux_basis,
                  const std::vector<double> &fitting,
                  const CPackedMatrix       &density,
                  const double               coulomb_factor,
                  const double               metric_factor,
                  const std::vector<int>    &atoms,
                  const std::vector<int>    &aux_atoms,
                  const bool                 with_metric) const -> CPackedMatrix;

    /// @brief The screening threshold of the integrals.
    double _threshold;

    /// @brief The target number of atom pairs of a block.
    size_t _block_size;
};

#endif /* SimdRIJGradientDriver_hpp */
