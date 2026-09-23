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
#include <utility>
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

/// @brief The fitted densities of a gradient divided over the ranks of a
/// communicator, which each rank holds a share of.
///
/// @note **The rows are the share and the columns are not.** A rank owns the
/// auxiliary functions of the atoms whose B vectors it holds, and holds the fitted
/// densities of those and of no others. `densities` is the length of the whole
/// auxiliary basis all the same, with only the owned entries allocated: an entry
/// nobody filled is an empty matrix and costs the header, so the contraction which
/// indexes it by the global auxiliary function needs no second index and no change.
///
/// @note The Gram rows are this rank's rows of a matrix the whole of which is
/// needed, and are gathered once at the end. It is the square of the auxiliary
/// basis, an order of magnitude below the fitted densities on anything worth
/// dividing, which is why it is the one thing here that ends up on every rank.
struct TDistributedFit
{
    /// @brief The auxiliary functions this rank owns, as indices into the whole
    /// auxiliary basis, ascending.
    std::vector<size_t> functions;

    /// @brief The fitted densities, the length of the auxiliary basis, with only
    /// the owned entries allocated.
    std::vector<CPackedMatrix> densities;

    /// @brief The fitted densities of the second spin, empty for a closed shell.
    std::vector<CPackedMatrix> densities_beta;

    /// @brief The fitting coefficients. Partial until the ranks have added theirs
    /// and the transposed factor has been applied to the sum.
    std::vector<double> fitting;

    /// @brief This rank's rows of the Gram product, its own functions by the whole
    /// auxiliary basis, accumulated one panel at a time.
    std::vector<double> gram_rows;

    /// @brief The dimensions of the auxiliary basis.
    size_t naux = 0;

    /// @brief The elements of one fitted density, the packed orbital pairs.
    /// @note **One of these for each spin.** The two spins of an open shell occupy
    /// different numbers of orbitals, so their fitted densities are different sizes
    /// and their panels are different panels. Running the second spin over the
    /// first's geometry reads past the end of every matrix it touches.
    size_t nelements = 0;

    /// @brief The elements a panel carries.
    size_t npanel = 0;

    /// @brief The elements of one fitted density of the second spin.
    size_t nelements_beta = 0;

    /// @brief The elements a panel of the second spin carries.
    size_t npanel_beta = 0;

    /// @brief The number of panels the elements of a spin are divided into.
    auto panels(const bool beta = false) const -> size_t
    {
        const auto total = beta ? nelements_beta : nelements;

        const auto width = beta ? npanel_beta : npanel;

        return (width == 0) ? 0 : (total + width - 1) / width;
    }

    /// @brief The first element of a panel of a spin and how many it carries.
    auto panel_range(const size_t ipanel, const bool beta = false) const -> std::pair<size_t, size_t>
    {
        const auto total = beta ? nelements_beta : nelements;

        const auto width = beta ? npanel_beta : npanel;

        const auto first = ipanel * width;

        return {first, std::min(width, (first < total) ? total - first : size_t{0})};
    }
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

    /// @brief Forms this rank's share of the fitted densities, with no metric
    /// applied and nothing communicated.
    /// @param bq_vectors The B vectors this rank holds, of its auxiliary atoms.
    /// @param coefficients The occupied molecular orbitals.
    /// @param aux_atoms The atoms of the auxiliary basis this rank holds.
    /// @param budget The memory this phase may hold, which under a communicator is
    /// a share of the node's and not the whole of it.
    /// @return The share, with `fitting` holding this rank's **partial** sum of the
    /// right hand side of the fitting, which the caller adds across the ranks.
    /// @note This is the expensive phase and the one which needed no arranging: the
    /// half transformed W vectors come from the B vectors of this rank's own atoms,
    /// so the transformation is divided already by the division the Fock build made.
    /// @param coefficients_beta The occupied orbitals of the second spin, or an
    /// empty matrix for a closed shell. An open shell fits a set of densities for
    /// each spin, from the same B vectors and its own orbitals.
    auto mpi_local_densities(const CMolecule        &molecule,
                             const CMolecularBasis  &basis,
                             const CMolecularBasis  &aux_basis,
                             const CSparseTensor    &bq_vectors,
                             const CPackedMatrix    &density,
                             const CPackedMatrix    &coefficients,
                             const CPackedMatrix    &coefficients_beta,
                             const std::vector<int> &aux_atoms,
                             const size_t            budget) const -> TDistributedFit;

    /// @brief Applies the transposed factor to the fitting coefficients the ranks
    /// have added together.
    /// @param fit The share, whose `fitting` is replaced by the coefficients.
    /// @param metric The inverted factor of the metric.
    /// @param total The sum over the ranks of the right hand side.
    auto mpi_set_fitting(TDistributedFit           &fit,
                         const CPackedMatrix       &metric,
                         const std::vector<double> &total) const -> void;

    /// @brief Empties the fitting coefficients of a share, which is how the
    /// attenuated operator says it has no Coulomb term.
    /// @note Empty and not zero, so that anything which reaches for them fails
    /// rather than returning a number which looks like an answer. It is the
    /// convention the serial attenuated phase already uses.
    auto mpi_clear_fitting(TDistributedFit &fit) const -> void;

    /// @brief This rank's contribution to one panel of the transposed factor
    /// applied to the fitted densities.
    /// @return The whole auxiliary basis by the elements of the panel, row major,
    /// holding the part of the sum over the auxiliary functions this rank owns.
    /// @note Every rank answers for every auxiliary function and for the elements
    /// of this panel alone. What makes the phase divisible is that the sum over the
    /// **owned** index is the one which factorises; the sum over the answered index
    /// is what the ranks then add.
    auto mpi_panel_partial(const TDistributedFit &fit,
                           const CPackedMatrix   &metric,
                           const size_t           ipanel,
                           const bool             beta = false) const -> std::vector<double>;

    /// @brief Takes one panel the ranks have added together, keeping the rows this
    /// rank owns and accumulating its rows of the Gram product from the whole.
    /// @param fit The share.
    /// @param ipanel The panel.
    /// @param reduced The panel, summed over the ranks.
    /// @param beta Whether the panel is the second spin's.
    /// @note The Gram is accumulated here rather than afterwards because this is
    /// the one moment the whole of a panel exists on a rank. Afterwards only the
    /// owned rows remain, and the Gram of a pair of rows on different ranks could
    /// not be formed at all.
    /// @note The panel is taken as a pointer and not as a vector, so that a caller
    /// handing over a buffer of the size the budget allows does not copy it first.
    auto mpi_panel_absorb(TDistributedFit &fit,
                          const size_t     ipanel,
                          const double    *reduced,
                          const size_t     size,
                          const bool       beta = false) const -> void;

    /// @brief Assembles the two-index fitted density from the fitting coefficients
    /// and the Gram rows the ranks have gathered.
    /// @param fit The share.
    /// @param gram The Gram, the whole of it, gathered from the rows of the ranks.
    /// @param exchange_scaling_factor The fraction of exact exchange.
    /// @param open_shell Whether the factors are an open shell's.
    /// @return Omega, of the whole auxiliary basis.
    auto mpi_omega(const TDistributedFit &fit,
                   const double          *gram,
                   const size_t           size,
                   const double           exchange_scaling_factor,
                   const bool             open_shell) const -> CPackedMatrix;

    /// @brief This rank's share of the gradient, from the fitted densities it owns.
    /// @param fit The share, with its densities transformed and its fitting set.
    /// @param omega The two-index fitted density, or an empty matrix to leave the
    /// two-center term out, which every rank but one does.
    /// @param aux_atoms The atoms of the auxiliary basis this rank holds.
    /// @return The gradient of this rank's share, which the ranks add.
    /// @param coefficients_beta The occupied orbitals of the second spin, empty
    /// for a closed shell. A spin occupies its own orbitals, and handing the first
    /// spin's for both is an exchange of a wavefunction nobody asked for.
    auto mpi_compute_share(const CMolecule        &molecule,
                           const CMolecularBasis  &basis,
                           const CMolecularBasis  &aux_basis,
                           const TDistributedFit  &fit,
                           const CPackedMatrix    &density,
                           const CPackedMatrix    &coefficients,
                           const CPackedMatrix    &coefficients_beta,
                           const CPackedMatrix    &omega,
                           const double            exchange_scaling_factor,
                           const std::vector<int> &atoms,
                           const std::vector<int> &aux_atoms) const -> CPackedMatrix;

    /// @brief This rank's share of a range separated gradient, from the plain and
    /// the attenuated fitted densities it owns.
    /// @note The attenuated operator has no Coulomb term, so its fit comes with an
    /// empty fitting and its Omega carries the exchange alone.
    auto mpi_compute_share_rs(const CMolecule        &molecule,
                              const CMolecularBasis  &basis,
                              const CMolecularBasis  &aux_basis,
                              const TDistributedFit  &fit,
                              const TDistributedFit  &fit_erf,
                              const CPackedMatrix    &density,
                              const CPackedMatrix    &coefficients,
                              const CPackedMatrix    &coefficients_beta,
                              const CPackedMatrix    &omega_plain,
                              const CPackedMatrix    &omega_erf,
                              const double            exchange_scaling_factor,
                              const double            erf_exchange_scaling_factor,
                              const double            omega,
                              const std::vector<int> &atoms,
                              const std::vector<int> &aux_atoms) const -> CPackedMatrix;

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
    /// @brief Computes the gradient of a hybrid range separated functional, the
    /// restricted form.
    /// @param bq_vectors The B vectors of the Coulomb operator.
    /// @param bq_vectors_erf The B vectors of the attenuated operator.
    /// @param metric The inverted metric of the Coulomb operator.
    /// @param metric_erf The inverted metric of the attenuated operator, which is
    /// a different matrix and not interchangeable with it.
    /// @param exchange_scaling_factor The coefficient of the plain exchange, which
    /// is alpha plus beta.
    /// @param erf_exchange_scaling_factor The coefficient of the attenuated
    /// exchange, which is the erf coefficient of the functional.
    /// @param omega The range separation parameter, which must be positive.
    /// @note A separate entry rather than a flag on compute, so that a calculation
    /// which is not range separated cannot reach this code and a reader of either
    /// can see which case they are in. The Coulomb term appears once and is fitted
    /// in the plain metric: only the exchange is split.
    auto compute_rs(const CMolecule        &molecule,
                    const CMolecularBasis  &basis,
                    const CMolecularBasis  &aux_basis,
                    const CSparseTensor    &bq_vectors,
                    const CSparseTensor    &bq_vectors_erf,
                    const CPackedMatrix    &metric,
                    const CPackedMatrix    &metric_erf,
                    const CPackedMatrix    &density,
                    const CPackedMatrix    &coefficients,
                    const double            exchange_scaling_factor,
                    const double            erf_exchange_scaling_factor,
                    const double            omega,
                    const std::vector<int> &atoms,
                    const std::vector<int> &aux_atoms = {}) const -> CPackedMatrix;

    /// @brief Computes it for every atom of the molecule.
    auto compute_rs(const CMolecule       &molecule,
                    const CMolecularBasis &basis,
                    const CMolecularBasis &aux_basis,
                    const CSparseTensor   &bq_vectors,
                    const CSparseTensor   &bq_vectors_erf,
                    const CPackedMatrix   &metric,
                    const CPackedMatrix   &metric_erf,
                    const CPackedMatrix   &density,
                    const CPackedMatrix   &coefficients,
                    const double           exchange_scaling_factor,
                    const double           erf_exchange_scaling_factor,
                    const double           omega) const -> CPackedMatrix;

    /// @brief Computes the gradient of a hybrid range separated functional, the
    /// unrestricted form.
    /// @param density The total density, which is that of both spins added.
    /// @param coefficients_alpha The occupied orbitals of the alpha spin.
    /// @param coefficients_beta The same for the beta spin, which has a number of
    /// columns of its own.
    /// @note The factors are those of compute_open_shell: the Coulomb carries one
    /// where the closed shell carries four, and each spin's exchange carries its
    /// coefficient once where the closed shell carries it twice.
    auto compute_open_shell_rs(const CMolecule        &molecule,
                               const CMolecularBasis  &basis,
                               const CMolecularBasis  &aux_basis,
                               const CSparseTensor    &bq_vectors,
                               const CSparseTensor    &bq_vectors_erf,
                               const CPackedMatrix    &metric,
                               const CPackedMatrix    &metric_erf,
                               const CPackedMatrix    &density,
                               const CPackedMatrix    &coefficients_alpha,
                               const CPackedMatrix    &coefficients_beta,
                               const double            exchange_scaling_factor,
                               const double            erf_exchange_scaling_factor,
                               const double            omega,
                               const std::vector<int> &atoms,
                               const std::vector<int> &aux_atoms = {}) const -> CPackedMatrix;

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
                               const std::vector<int>           &aux_atoms,
                               const std::vector<TExchangeSpin> &spins_erf = {},
                               const double                      erf_exchange_factor = 0.0,
                               const double                      omega = 0.0) const -> void;

    /// @brief Forms one spin's fitted densities for the functions this rank owns.
    auto _local_densities_for(const CSparseTensor       &bq_vectors,
                              const CMolecularBasis     &basis,
                              const CMolecularBasis     &aux_basis,
                              const CPackedMatrix       &coefficients,
                              const std::vector<size_t> &functions,
                              const size_t               naux,
                              const size_t               budget,
                              std::vector<CPackedMatrix> &target) const -> void;

    /// @brief The columns of the transposed factor of the metric which belong to
    /// the given auxiliary functions.
    /// @param metric The inverted factor of the metric.
    /// @param columns The auxiliary functions whose columns are wanted.
    /// @param naux The dimensions of the auxiliary basis.
    /// @return The whole auxiliary basis by those columns, row major.
    auto _transposed_columns(const CPackedMatrix       &metric,
                             const std::vector<size_t> &columns,
                             const size_t               naux) const -> std::vector<double>;

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
