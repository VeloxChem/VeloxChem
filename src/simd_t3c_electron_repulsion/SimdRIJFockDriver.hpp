//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//

#ifndef SimdRIJFockDriver_hpp
#define SimdRIJFockDriver_hpp

#include <cstddef>
#include <vector>

#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "PackedMatrix.hpp"
#include "SimdRIFockCommon.hpp"
#include "SimdRIFockDriver.hpp"
#include "SparseTensor.hpp"
#include "TripleSparsityPattern.hpp"

/// @brief The Coulomb matrix of the resolution of the identity, for a functional
/// which asks for no exact exchange.
///
/// @note This is three steps and the middle one is why it is cheap:
///
///     g(P)     = sum over mn of (mn|P) D(mn)
///     gamma    = J inverse g
///     J(mn)    = sum over P of (mn|P) gamma(P)
///
/// **The metric never multiplies the tensor.** It multiplies a vector of one value
/// per auxiliary function, twice a build, which is the square of the auxiliary basis
/// against one right hand side and is nothing beside either sweep. The driver which
/// also forms the exchange cannot do this: an exchange needs the fitting closed for
/// every orbital, so it folds a factor of the metric into the integrals themselves
/// and pays the square of the auxiliary basis for every pair of orbital functions.
/// That transformation is ninety per cent of what its setup costs and none of it is
/// here.
///
/// @note So a pure functional wants this driver and a hybrid wants the other one,
/// and the difference between them is not the exchange it adds but the work it does
/// not have to do to get the Coulomb.
class CSimdRIJFockDriver
{
   public:
    /// @brief The default constructor.
    CSimdRIJFockDriver() = default;

    CSimdRIJFockDriver(const CSimdRIJFockDriver &) = delete;
    CSimdRIJFockDriver(CSimdRIJFockDriver &&) noexcept = delete;
    ~CSimdRIJFockDriver() = default;
    auto operator=(const CSimdRIJFockDriver &) -> CSimdRIJFockDriver & = delete;
    auto operator=(CSimdRIJFockDriver &&) noexcept -> CSimdRIJFockDriver & = delete;

    /// @brief Gets the memory the driver would hold for a molecule and its bases.
    /// @param aux_atoms The atoms of the auxiliary side this rank holds, or none of
    /// them for the whole molecule.
    /// @return The memory, in bytes.
    /// @note One tensor, where the driver which forms the exchange holds one or two
    /// of them and says so. Nothing else is held: there is no transformation of the
    /// integrals to keep.
    auto required_memory(const CMolecule        &molecule,
                         const CMolecularBasis  &basis,
                         const CMolecularBasis  &aux_basis,
                         const double            threshold,
                         const std::vector<int> &aux_atoms = {}) const -> size_t;

    /// @brief Forms the metric of the auxiliary basis, inverted outright.
    /// @param metric_threshold The threshold of the linear dependence.
    /// @return The inverted metric.
    /// @note Formed once by whoever prepares the calculation and handed to prepare,
    /// so that the ranks of a communicator build with the same matrix rather than
    /// each inverting and risking a different fallback.
    auto make_metric(const CMolecule       &molecule,
                     const CMolecularBasis &aux_basis,
                     const double           metric_threshold) const -> CPackedMatrix;

    /// @brief Prepares the driver for a molecule and its bases.
    /// @param threshold The screening threshold of the integrals.
    /// @param memory_budget The memory a part's integrals may take, in bytes.
    /// @param mode The way to build. Automatic chooses from the memory.
    /// @param metric The inverted metric, or an empty matrix to form it here.
    /// @param rank This rank, which decides which parts it owns.
    /// @param nodes The ranks which divide the parts, so that there are at least as
    /// many parts as there are of them.
    /// @note **A rank holds only the parts it sweeps.** The two sweeps are divided
    /// over the parts, so a rank never reads another's; holding them all would cost
    /// the memory of the whole molecule on every one of them, and would send the
    /// automatic choice to the direct way for a calculation each rank holds a
    /// fitting share of.
    auto prepare(const CMolecule       &molecule,
                 const CMolecularBasis &basis,
                 const CMolecularBasis &aux_basis,
                 const double           threshold,
                 const size_t           memory_budget,
                 const double           metric_threshold,
                 const rimode           mode,
                 const CPackedMatrix   &metric,
                 const size_t           rank,
                 const size_t           nodes) -> void;

    /// @brief Checks that the driver has been prepared.
    auto is_prepared() const -> bool;

    /// @brief Gets the way the driver builds, which is never automatic once prepared.
    auto get_mode() const -> rimode;

    /// @brief Gets the inverted metric the driver holds.
    auto get_metric() const -> const CPackedMatrix &;

    /// @brief Gets the parts of the auxiliary basis the sweeps take.
    /// @return The number of parts, which bounds the indices the sweeps accept.
    auto number_of_parts() const -> size_t;

    /// @brief Gets the parts this rank owns, which are the ones it is to sweep.
    /// @return The indices, as the sweeps take them.
    /// @note Asked of the driver rather than worked out again by the caller. The
    /// driver held these parts and no others, so a caller which divided them its own
    /// way would ask a rank for a part whose integrals it does not have.
    auto owned_parts() const -> std::vector<int>;

    /// @brief Computes the Coulomb matrix of a density, on one rank.
    /// @param density The density matrix.
    /// @return The Coulomb matrix, J(mn) = sum over P of (mn|P) gamma(P).
    /// @note **Once and not twice.** The caller of a restricted calculation hands
    /// the density of one spin and doubles what comes back, which is what the
    /// conventional driver of this approximation does and what the Fock build
    /// expects of a matrix of the `j` kind.
    auto compute(const CPackedMatrix &density) -> CPackedMatrix;

    /// @brief Sweeps the given parts and sums the right hand side of the fitting.
    /// @param parts The parts to sweep, as their indices.
    /// @return The right hand side, of one value per auxiliary basis function, of
    /// those parts alone.
    /// @note The ranks of a communicator take some parts each and **add theirs
    /// together** before the fitting: the solve reaches across the whole auxiliary
    /// basis and needs a right hand side which is complete.
    auto compute_gamma(const CPackedMatrix &density, const std::vector<int> &parts) -> std::vector<double>;

    /// @brief Applies the inverted metric to the right hand side of the fitting.
    /// @param gamma The right hand side, summed over every part there is.
    /// @return The coefficients of the fitting.
    /// @note Every rank holds the inverse and solves it, rather than one solving and
    /// sending: it is one multiply of the square of the auxiliary basis against a
    /// single right hand side.
    auto solve_fitting(const std::vector<double> &gamma) const -> std::vector<double>;

    /// @brief Adds the Coulomb matrix of the given parts to a matrix.
    /// @param gamma The coefficients of the fitting.
    /// @param parts The parts to sweep, as their indices.
    /// @param matrix The matrix to add to.
    auto compute_coulomb(const std::vector<double> &gamma, const std::vector<int> &parts, CPackedMatrix &matrix) -> void;

   private:
    /// @brief Points at the integrals of one part, forming them if they are not held.
    /// @param index The part.
    /// @param formed Where to form them when they are not held, which the caller owns
    /// and which is left alone when they are.
    /// @return A pointer to the integrals, either the held ones or the formed ones.
    /// @note A pointer and not a value. Handing back a copy would give the two ways
    /// one signature, and it would copy the whole of a part on every sweep of every
    /// iteration: at 32 waters that was most of the gap between this driver and the
    /// conventional one, which holds its integrals and reads them where they are.
    auto _part_integrals(const size_t index, CSparseTensor &formed) -> const CSparseTensor *;

    /// @brief Checks that an index names a part this driver sweeps.
    auto _check_part(const int index) const -> void;

    bool _prepared = false;

    rimode _mode = rimode::automatic;

    CMolecule _molecule;

    CMolecularBasis _basis;

    CMolecularBasis _aux_basis;

    /// @brief The metric of the auxiliary basis, inverted outright.
    CPackedMatrix _metric;

    /// @brief The sparsity pattern of each part of the auxiliary basis.
    std::vector<CTripleSparsityPattern> _parts;

    /// @brief The integrals of each part, held only by the way which holds them and
    /// only for the parts this rank owns. The entries of the others stay empty, so
    /// that an index means the same thing here as it does in the parts.
    std::vector<CSparseTensor> _integrals;

    /// @brief The parts this rank owns.
    std::vector<int> _owned;

    /// @brief The driver of the sweeps, which the two ways share.
    CSimdRIFockDriver _drv;
};

#endif /* SimdRIJFockDriver_hpp */
