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

#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "PackedMatrix.hpp"
#include "SimdRIFockDriver.hpp"
#include "SimdThreeCenterElectronRepulsionGradientDriver.hpp"
#include "SparseTensor.hpp"

/// @brief The fitted densities the gradient contracts the derivative integrals
/// against, which are formed once for the whole gradient.
struct TFittedDensities
{
    /// @brief The fitting coefficients, c = V inverted times gamma, one per
    /// auxiliary basis function.
    std::vector<double> coefficients;

    /// @brief The fitted densities of the occupied orbitals, d(q)_ij, one packed
    /// symmetric matrix of the occupied orbitals for each auxiliary function.
    std::vector<CPackedMatrix> orbital_densities;

    /// @brief The fitted densities of the second spin, empty for a closed shell.
    std::vector<CPackedMatrix> orbital_densities_beta;

    /// @brief The two-index fitted density, one row and column per auxiliary
    /// basis function.
    CPackedMatrix omega;
};

/// @brief One spin's half of the exchange: the orbitals it occupies and the
/// fitted densities formed from them.
/// @note A closed shell is one of these and an open shell two. The exchange of an
/// auxiliary function is their sum, which is why they travel together.
struct TExchangeSpin
{
    const CPackedMatrix              *coefficients;
    const std::vector<CPackedMatrix> *orbital_densities;
};

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
    /// @param memory_budget The memory the driver may hold for the half
    /// transformed B vectors, in bytes. It bounds the batch of auxiliary
    /// functions and nothing else: the fitted densities of the orbitals are
    /// M o squared over two whatever the batch, and are held whole.
    explicit CSimdRIJKGradientDriver(const double threshold      = 1.0e-12,
                                     const size_t block_size     = 0,
                                     const size_t memory_budget  = _default_budget)
        : _threshold(threshold)
        , _block_size(block_size)
        , _budget(memory_budget)
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

    /// @brief Gets the memory the driver may hold, in bytes.
    /// @return The memory budget.
    auto get_memory_budget() const -> size_t;

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
    /// @param exchange_scaling_factor The fraction of exact exchange: one for
    /// Hartree-Fock, the fraction of a hybrid, zero for a pure functional.
    /// @param atoms The atoms to compute the gradient of.
    /// @param aux_atoms The atoms of the auxiliary basis the B vectors span, or
    /// an empty list for all of them.
    /// @return The gradient, a general matrix of one row of three components per
    /// atom of the molecule, with the rows of the atoms not asked for left zero.
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
                 const double            exchange_scaling_factor,
                 const std::vector<int> &atoms,
                 const std::vector<int> &aux_atoms = {}) const -> CPackedMatrix;

    /// @brief Computes the Coulomb and exchange contributions to the gradient of
    /// every atom of the molecule.
    /// @return The gradient, a general matrix of three components per atom.
    auto compute(const CMolecule       &molecule,
                 const CMolecularBasis &basis,
                 const CMolecularBasis &aux_basis,
                 const CSparseTensor   &bq_vectors,
                 const CPackedMatrix   &metric,
                 const CPackedMatrix   &density,
                 const CPackedMatrix   &coefficients,
                 const double           exchange_scaling_factor) const -> CPackedMatrix;

    /// @brief The fitted densities of the occupied orbitals of one spin.
    /// @param bq_vectors The B vectors.
    /// @param basis The molecular basis.
    /// @param aux_basis The auxiliary molecular basis.
    /// @param metric The inverted Cholesky factor of the metric.
    /// @param coefficients The occupied orbitals of the spin.
    /// @return d(q)_ij with the metric applied, one packed matrix per auxiliary
    /// function.
    /// @note One spin's worth, so that an open shell forms two sets from one
    /// routine rather than two copies of a loop which must agree.
    auto _orbital_densities(const CSparseTensor   &bq_vectors,
                            const CMolecularBasis &basis,
                            const CMolecularBasis &aux_basis,
                            const CPackedMatrix   &metric,
                            const CPackedMatrix   &coefficients) const -> std::vector<CPackedMatrix>;

    /// @brief The Gram product of the fitted densities over the pairs of
    /// auxiliary functions, summed over the orbitals.
    /// @param orbital_densities The fitted densities of one spin.
    /// @param naux The number of auxiliary functions.
    /// @param norbs The number of occupied orbitals of that spin.
    /// @return The square of the auxiliary basis, row major.
    auto _gram(const std::vector<CPackedMatrix> &orbital_densities,
               const size_t                      naux,
               const size_t                      norbs) const -> std::vector<double>;

    /// @brief Computes the gradient of an open shell, the two spins occupying
    /// different orbitals.
    /// @param density The **total** density, of both spins added, which is what
    /// the Coulomb half is of. The closed shell call above takes one spin's and
    /// carries the factors of two which follow from that.
    /// @param coefficients_alpha The occupied orbitals of the alpha spin.
    /// @param coefficients_beta The occupied orbitals of the beta spin, which has
    /// a number of its own.
    /// @param exchange_scaling_factor The fraction of exact exchange.
    /// @param atoms The atoms to compute the gradient of.
    /// @param aux_atoms The atoms of the auxiliary basis the B vectors span.
    /// @return The gradient, one row per atom.
    /// @note A separate routine and not an overload of the one above: the two
    /// differ in which density the Coulomb half is of, and a caller which passed
    /// the wrong one would get a gradient that is merely wrong. The name says
    /// which is meant.
    auto compute_open_shell(const CMolecule        &molecule,
                            const CMolecularBasis  &basis,
                            const CMolecularBasis  &aux_basis,
                            const CSparseTensor    &bq_vectors,
                            const CPackedMatrix    &metric,
                            const CPackedMatrix    &density,
                            const CPackedMatrix    &coefficients_alpha,
                            const CPackedMatrix    &coefficients_beta,
                            const double            exchange_scaling_factor,
                            const std::vector<int> &atoms,
                            const std::vector<int> &aux_atoms) const -> CPackedMatrix;

    /// @brief Computes the gradient of an open shell over every atom.
    auto compute_open_shell(const CMolecule       &molecule,
                            const CMolecularBasis &basis,
                            const CMolecularBasis &aux_basis,
                            const CSparseTensor   &bq_vectors,
                            const CPackedMatrix   &metric,
                            const CPackedMatrix   &density,
                            const CPackedMatrix   &coefficients_alpha,
                            const CPackedMatrix   &coefficients_beta,
                            const double           exchange_scaling_factor) const -> CPackedMatrix;

    /// @brief Forms the fitted densities of an open shell, both spins.
    /// @param density The total density, which the fitting coefficients are of.
    /// @param coefficients_alpha The occupied orbitals of the alpha spin.
    /// @param coefficients_beta The occupied orbitals of the beta spin.
    /// @param exchange_scaling_factor The fraction of exact exchange.
    /// @return The fitting coefficients, a set of fitted densities for each spin,
    /// and the two index fitted density of both.
    auto fitted_densities_open_shell(const CSparseTensor   &bq_vectors,
                                     const CMolecularBasis &basis,
                                     const CMolecularBasis &aux_basis,
                                     const CPackedMatrix   &metric,
                                     const CPackedMatrix   &density,
                                     const CPackedMatrix   &coefficients_alpha,
                                     const CPackedMatrix   &coefficients_beta,
                                     const double           exchange_scaling_factor) const -> TFittedDensities;

    /// @brief Forms the fitted densities the derivative integrals are contracted
    /// against, which is the whole of what the gradient needs before any of them
    /// are computed.
    /// @param bq_vectors The B vectors, as CSimdRIJKFockDriver formed them.
    /// @param basis The molecular basis.
    /// @param aux_basis The auxiliary molecular basis.
    /// @param metric The inverted Cholesky factor of the metric.
    /// @param density The density matrix.
    /// @param coefficients The occupied molecular orbitals.
    /// @param exchange_scaling_factor The fraction of exact exchange.
    /// @return The fitting coefficients, the fitted densities of the orbitals and
    /// the two-index fitted density.
    /// @note This phase reads integrals of none: every quantity here is a
    /// transformation of the B vectors the calculation already holds. The
    /// transformation into the occupied orbitals is taken **before** the metric is
    /// applied, which is the ordering rule of the note: applying the metric in the
    /// basis of the atomic orbitals costs the ratio of their number squared to the
    /// orbitals squared, which is a factor of fifty for a molecule of this size.
    /// @brief The fitted densities of the attenuated operator, for a hybrid range
    /// separated functional.
    /// @param bq_vectors_erf The B vectors of the attenuated operator.
    /// @param basis The molecular basis.
    /// @param aux_basis The auxiliary molecular basis.
    /// @param metric_erf The inverted factor of the attenuated metric. The two
    /// operators are fitted in metrics of their own and crossing them fits
    /// neither.
    /// @param coefficients The occupied molecular orbitals.
    /// @param erf_exchange_scaling_factor The coefficient of the attenuated
    /// exchange, which is the erf coefficient of the functional.
    /// @return The fitted densities of the orbitals and the two-index fitted
    /// density of the exchange. **The fitting coefficients come back empty**: the
    /// attenuated operator enters the Fock matrix only through the exchange, so it
    /// has no Coulomb term and nothing to fit the density against. A caller which
    /// hands these to the three-center term must pass a Coulomb factor of zero.
    auto fitted_densities_rs(const CSparseTensor   &bq_vectors_erf,
                             const CMolecularBasis &basis,
                             const CMolecularBasis &aux_basis,
                             const CPackedMatrix   &metric_erf,
                             const CPackedMatrix   &coefficients,
                             const double           erf_exchange_scaling_factor) const -> TFittedDensities;

    /// @brief The same for the two spins of an open shell.
    /// @param coefficients_alpha The occupied orbitals of the alpha spin.
    /// @param coefficients_beta The same for the beta spin, which has a number of
    /// columns of its own.
    /// @note The exchange carries a half here as it does in the unattenuated open
    /// shell routine, so that setting the two spins equal returns the closed shell
    /// expression exactly.
    auto fitted_densities_open_shell_rs(const CSparseTensor   &bq_vectors_erf,
                                        const CMolecularBasis &basis,
                                        const CMolecularBasis &aux_basis,
                                        const CPackedMatrix   &metric_erf,
                                        const CPackedMatrix   &coefficients_alpha,
                                        const CPackedMatrix   &coefficients_beta,
                                        const double erf_exchange_scaling_factor) const -> TFittedDensities;

    auto fitted_densities(const CSparseTensor   &bq_vectors,
                          const CMolecularBasis &basis,
                          const CMolecularBasis &aux_basis,
                          const CPackedMatrix   &metric,
                          const CPackedMatrix   &density,
                          const CPackedMatrix   &coefficients,
                          const double           exchange_scaling_factor) const -> TFittedDensities;

   private:
    /// @brief Applies the transpose of the inverted Cholesky factor to columns.
    auto _apply_transposed_factor(const CPackedMatrix &metric, double *values, const size_t ncols) const -> void;

    /// @brief Applies it over the auxiliary index of a set of matrices.
    auto _apply_transposed_factor(const CPackedMatrix &metric, std::vector<CPackedMatrix> &matrices) const -> void;

    /// @brief Closes the second index of a half transformed B vector into the
    /// occupied orbitals.
    auto _close_orbitals(const std::vector<double> &transposed,
                         const size_t               nao,
                         const size_t               norbs,
                         const CPackedMatrix       &half) const -> CPackedMatrix;

    /// @brief Adds the three-center term to the gradient, one atom of the
    /// auxiliary basis at a time.
    /// @note The derivative integrals of an auxiliary atom are formed, contracted
    /// against Gamma and discarded before the next atom is asked for. The Gamma
    /// of the atomic orbitals is never formed for the whole auxiliary basis: at
    /// the functions squared times the auxiliary basis it is hundreds of
    /// gigabytes, where one atom's share is the functions squared times the
    /// functions of that atom.
    auto _compute_three_center(CPackedMatrix                    &gradient,
                               const CMolecule                  &molecule,
                               const CMolecularBasis            &basis,
                               const CMolecularBasis            &aux_basis,
                               const std::vector<double>        &fitting,
                               const CPackedMatrix              &density,
                               const std::vector<TExchangeSpin> &spins,
                               const double                      coulomb_factor,
                               const double                      exchange_factor,
                               const std::vector<bool>          &wanted,
                               const std::vector<int>           &aux_atoms) const -> void;

    /// @brief Checks the metric is one this driver can use.
    /// @param metric The metric handed over.
    /// @param aux_basis The auxiliary molecular basis.
    auto _check_metric(const CPackedMatrix &metric, const CMolecularBasis &aux_basis) const -> void;

    /// @brief The driver which forms and transforms the B vectors.
    CSimdRIFockDriver _drv;

    /// @brief The screening threshold of the integrals.
    double _threshold;

    /// @brief The target number of atom pairs of a block.
    size_t _block_size;

    /// @brief The memory the driver may hold for the half transformed B vectors.
    size_t _budget;

    /// @brief The memory a driver may hold when none is named, in bytes.
    /// @note Four gigabytes, which is a batch of the auxiliary basis and not the
    /// fitted densities themselves. A caller which knows what the machine has, or
    /// how many of it share the machine, should say so rather than take this.
    static constexpr size_t _default_budget = size_t{4} * 1024 * 1024 * 1024;

    /// @brief The smallest batch of auxiliary functions to transform at a time.
    /// @note A batch below this is not worth the call, and one function at a time
    /// is what the test which checks the batching does not change the answer asks
    /// for, so the floor is not applied when the budget names something smaller.
    static constexpr size_t _min_batch = 16;
};

#endif /* SimdRIJKGradientDriver_hpp */
