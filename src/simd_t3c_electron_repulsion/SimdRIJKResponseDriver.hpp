//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//

#ifndef SimdRIJKResponseDriver_hpp
#define SimdRIJKResponseDriver_hpp

#include <cstddef>
#include <utility>
#include <vector>

#include "MolecularBasis.hpp"
#include "PackedMatrix.hpp"
#include "SimdRIFockDriver.hpp"
#include "SparseTensor.hpp"

/// @brief The Fock matrices of a response calculation, through the resolution of
/// the identity, from B vectors which are already formed.
///
/// @note This driver is not the one the SCF uses and shares none of its routines.
/// The two are asked for different things. The SCF has one density per build, it
/// is symmetric, and it is the product of the occupied orbitals with themselves,
/// which is what lets its exchange be formed from the orbitals alone. A response
/// calculation has many densities per build, none of them symmetric and none of
/// them idempotent: a trial vector gives a density which lives entirely between
/// the occupied orbitals and the virtual ones. Neither the argument list nor the
/// exchange of the SCF driver applies, so they are kept apart rather than one of
/// them made to answer both.
///
/// @note The densities are taken in the factorised form they are made in and not
/// as matrices. A trial vector of the Tamm-Dancoff approximation gives the density
/// C(occupied) Z C(virtual) transposed, and the driver is handed those factors
/// rather than their product: the exchange of a density of rank k costs the
/// auxiliary basis times the square of the basis times k, where the exchange of
/// the same density as a matrix costs the auxiliary basis times the cube of the
/// basis. Forming the product and factorising it again would be the more
/// expensive of the two, and the caller has the factors already.
class CSimdRIJKResponseDriver
{
   public:
    /// @brief Creates a response driver with a screening threshold and a target
    /// block size.
    /// @param threshold The screening threshold of the integrals.
    /// @param block_size The target number of atom pairs of a block, or zero to
    /// choose it from the number of the threads and the number of the atom pairs.
    /// @param memory_budget The memory the driver may hold for the transformed B
    /// vectors, in bytes. It bounds the batch of auxiliary functions and nothing
    /// else.
    explicit CSimdRIJKResponseDriver(const double threshold     = 1.0e-12,
                                     const size_t block_size    = 0,
                                     const size_t memory_budget = _default_budget)
        : _threshold(threshold)
        , _block_size(block_size)
        , _budget(memory_budget)
    {
    }

    CSimdRIJKResponseDriver(const CSimdRIJKResponseDriver &other)                = delete;
    CSimdRIJKResponseDriver(CSimdRIJKResponseDriver &&other) noexcept            = delete;
    ~CSimdRIJKResponseDriver()                                                   = default;
    auto operator=(const CSimdRIJKResponseDriver &other) -> CSimdRIJKResponseDriver     & = delete;
    auto operator=(CSimdRIJKResponseDriver &&other) noexcept -> CSimdRIJKResponseDriver & = delete;

    /// @brief Gets the screening threshold of the integrals.
    /// @return The screening threshold.
    auto get_threshold() const -> double;

    /// @brief Gets the target number of atom pairs of a block.
    /// @return The target number of atom pairs.
    auto get_block_size() const -> size_t;

    /// @brief Gets the memory the driver may hold, in bytes.
    /// @return The memory budget.
    auto get_memory_budget() const -> size_t;

    /// @brief Computes the exchange matrices of densities given as their factors.
    /// @param bq_vectors The B vectors, as CSimdRIJKFockDriver formed them.
    /// @param basis The molecular basis.
    /// @param aux_basis The auxiliary molecular basis the identity is resolved in.
    /// @param left The left factor, which every density of the batch shares.
    /// @param rights The right factors, one for each density of the batch. The
    /// density of one of them is the left factor times its transpose.
    /// @return One exchange matrix for each right factor, general and not
    /// symmetric, in the order the right factors were given.
    /// @note The sum which is formed is, for each density,
    /// K = sum over q of (B(q) left) (B(q) right) transposed, which is the
    /// exchange of the density left times right transposed. The left factor is
    /// transformed once for the whole batch: it is the occupied orbitals of the
    /// ground state for every trial vector a response calculation makes, and the
    /// batch is where that is worth exploiting.
    /// @note What is returned carries no scaling. The fraction of exact exchange
    /// belongs to the caller, which knows what it is asking for.
    auto compute_exchange(const CSparseTensor              &bq_vectors,
                          const CMolecularBasis            &basis,
                          const CMolecularBasis            &aux_basis,
                          const CPackedMatrix              &left,
                          const std::vector<CPackedMatrix> &rights) const -> std::vector<CPackedMatrix>;

    /// @brief Computes the exchange of the Coulomb operator and of the attenuated
    /// one together, each scaled, for a hybrid range separated functional.
    /// @param bq_vectors The B vectors of the Coulomb operator.
    /// @param bq_vectors_erf The B vectors of the attenuated operator, on the same
    /// sparsity pattern, as CSimdRIJKFockDriver formed them together.
    /// @param exchange_scaling_factor What the plain exchange is scaled by.
    /// @param erf_exchange_scaling_factor What the attenuated exchange is scaled by.
    /// @return One matrix for each right factor, holding the two exchanges already
    /// scaled and added. **Unlike the form above, what is returned carries its
    /// scaling**, because keeping the operators apart would mean a second matrix of
    /// the basis squared for every density of the batch, to be added and dropped.
    /// @note Both operators are swept inside one pass over the auxiliary basis and
    /// accumulate into one set of matrices. Each of them has its own B vectors and
    /// so its own transformation of the factors, which cannot be shared; the
    /// batching, the stacking of the right factors and the output are.
    auto compute_exchange(const CSparseTensor              &bq_vectors,
                          const CSparseTensor              &bq_vectors_erf,
                          const CMolecularBasis            &basis,
                          const CMolecularBasis            &aux_basis,
                          const CPackedMatrix              &left,
                          const std::vector<CPackedMatrix> &rights,
                          const double                      exchange_scaling_factor,
                          const double erf_exchange_scaling_factor) const -> std::vector<CPackedMatrix>;

    /// @brief Computes the Fock matrices of densities given as their factors.
    /// @param bq_vectors The B vectors, as CSimdRIJKFockDriver formed them.
    /// @param basis The molecular basis.
    /// @param aux_basis The auxiliary molecular basis the identity is resolved in.
    /// @param left The left factor, which every density of the batch shares.
    /// @param rights The right factors, one for each density of the batch.
    /// @param exchange_scaling_factor The fraction of exact exchange: one for
    /// Hartree-Fock, the fraction of a hybrid, zero for a pure functional.
    /// @return One Fock matrix for each right factor, twice the Coulomb less the
    /// scaled exchange, general and not symmetric.
    /// @note The matrix returned is what the response builders call 2jk and 2jkx,
    /// assembled here rather than by the caller. A pure functional asks for a
    /// scaling of zero and is given twice the Coulomb, which is what that path
    /// forms for itself by doubling what it is handed.
    /// @note The Coulomb closes the whole density and not its factors, so the
    /// density is formed from them here. It is the basis squared by the rank for
    /// each density, which is nothing beside the transformation of the B vectors,
    /// and the alternative is to ask the caller for a matrix it would have to
    /// build from the same two factors.
    auto compute(const CSparseTensor              &bq_vectors,
                 const CMolecularBasis            &basis,
                 const CMolecularBasis            &aux_basis,
                 const CPackedMatrix              &left,
                 const std::vector<CPackedMatrix> &rights,
                 const double                      exchange_scaling_factor,
                 const CSparseTensor              *bq_vectors_erf              = nullptr,
                 const double                      erf_exchange_scaling_factor = 0.0) const
        -> std::vector<CPackedMatrix>;

    /// @brief Computes them for densities of two terms, the second with the
    /// shared factor on the other side.
    /// @param rights The right factor of the first term of each density.
    /// @param transposed_rights The left factor of the second term of each
    /// density, which stands to the left of the shared factor transposed. Empty
    /// where the densities have one term, and then this is the form above.
    /// @note The density of one pair is left times rights transposed, plus
    /// transposed_rights times left transposed. That is the density a linear
    /// response trial vector makes: the excitation part carries the occupied
    /// orbitals on the left and the de-excitation part carries them on the right,
    /// and neither of them is the transpose of the other.
    /// @note Both terms are formed from one transformation of the shared factor.
    /// The exchange of the second is the exchange of its transpose transposed,
    /// which is why the virtual orbitals are never transformed; and the Coulomb
    /// sees only the symmetric part of a density, which the two terms share with
    /// the single term whose right factor is the sum of the two, so it is taken
    /// once and not twice.
    auto compute(const CSparseTensor              &bq_vectors,
                 const CMolecularBasis            &basis,
                 const CMolecularBasis            &aux_basis,
                 const CPackedMatrix              &left,
                 const std::vector<CPackedMatrix> &rights,
                 const std::vector<CPackedMatrix> &transposed_rights,
                 const double                      exchange_scaling_factor,
                 const CSparseTensor              *bq_vectors_erf              = nullptr,
                 const double                      erf_exchange_scaling_factor = 0.0) const
        -> std::vector<CPackedMatrix>;

    /// @brief Computes the Fock matrices of the two spins of an unrestricted
    /// reference, from densities given as the factors of each spin.
    /// @param bq_vectors The B vectors.
    /// @param basis The molecular basis.
    /// @param aux_basis The auxiliary molecular basis.
    /// @param left_alpha The left factor the alpha densities share, which is the
    /// occupied orbitals of that spin.
    /// @param rights_alpha The right factor of the first term of each alpha
    /// density.
    /// @param transposed_rights_alpha The left factor of the second term of each
    /// alpha density, empty where the densities have one term.
    /// @param left_beta The same for the beta spin, which has a number of columns
    /// of its own: the two spins of an open shell do not occupy the same number of
    /// orbitals.
    /// @param rights_beta As above, for the beta spin.
    /// @param transposed_rights_beta As above, for the beta spin.
    /// @param exchange_scaling_factor The fraction of exact exchange.
    /// @return The Fock matrices of the alpha spin and those of the beta spin, one
    /// of each per density.
    /// @note The Coulomb matrix is of the two spins' densities added and is **not**
    /// doubled, where the restricted entries above are handed one spin's density
    /// and double it. It is formed once for a pair and serves both, which is what
    /// makes this cheaper than two restricted builds.
    /// @note The exchange of each spin is formed from its own left factor. The two
    /// spins share no transformation: their left factors are the occupied orbitals
    /// of each of them and differ both in what they are and in how many there are.
    /// That is what an unrestricted batch costs over a restricted one.
    auto compute_unrestricted(const CSparseTensor              &bq_vectors,
                              const CMolecularBasis            &basis,
                              const CMolecularBasis            &aux_basis,
                              const CPackedMatrix              &left_alpha,
                              const std::vector<CPackedMatrix> &rights_alpha,
                              const std::vector<CPackedMatrix> &transposed_rights_alpha,
                              const CPackedMatrix              &left_beta,
                              const std::vector<CPackedMatrix> &rights_beta,
                              const std::vector<CPackedMatrix> &transposed_rights_beta,
                              const double                      exchange_scaling_factor) const
        -> std::pair<std::vector<CPackedMatrix>, std::vector<CPackedMatrix>>;

    /// @brief The same for a hybrid range separated functional, whose exchange is
    /// split between the plain operator and the attenuated one.
    /// @param bq_vectors_erf The B vectors of the attenuated operator, on the same
    /// sparsity pattern as the plain ones.
    /// @param erf_exchange_scaling_factor What the attenuated exchange is scaled
    /// by, which is the erf coefficient of the functional.
    /// @note A separate entry rather than a defaulted argument on the one above,
    /// so that a calculation which is not range separated cannot reach this code at
    /// all and a reader of either can see which case it is in.
    /// @note The Coulomb is formed from the plain B vectors here as everywhere: the
    /// splitting is of the exchange alone.
    auto compute_unrestricted_rs(const CSparseTensor              &bq_vectors,
                                 const CSparseTensor              &bq_vectors_erf,
                                 const CMolecularBasis            &basis,
                                 const CMolecularBasis            &aux_basis,
                                 const CPackedMatrix              &left_alpha,
                                 const std::vector<CPackedMatrix> &rights_alpha,
                                 const std::vector<CPackedMatrix> &transposed_rights_alpha,
                                 const CPackedMatrix              &left_beta,
                                 const std::vector<CPackedMatrix> &rights_beta,
                                 const std::vector<CPackedMatrix> &transposed_rights_beta,
                                 const double                      exchange_scaling_factor,
                                 const double                      erf_exchange_scaling_factor) const
        -> std::pair<std::vector<CPackedMatrix>, std::vector<CPackedMatrix>>;

   private:
    /// @brief Computes the exchange of one operator or of several, each scaled,
    /// into one set of matrices.
    /// @param operators The B vectors of each operator and what its exchange is
    /// scaled by. An operator whose scaling is zero is not swept at all.
    /// @param basis The molecular basis.
    /// @param aux_basis The auxiliary molecular basis.
    /// @param left The left factor, which every density of the batch shares.
    /// @param rights The right factors, one for each density of the batch.
    /// @return One matrix for each right factor.
    /// @note Both public forms above are written in terms of this, so the one
    /// operator case is the same code as the two operator one rather than a copy of
    /// it which drifts.
    auto _exchange(const std::vector<std::pair<const CSparseTensor *, double>> &operators,
                   const CMolecularBasis                                      &basis,
                   const CMolecularBasis                                      &aux_basis,
                   const CPackedMatrix                                        &left,
                   const std::vector<CPackedMatrix> &rights) const -> std::vector<CPackedMatrix>;

    /// @brief Computes the Fock matrices of the two spins, for one operator or
    /// several.
    /// @param operators The B vectors of each operator of the exchange and what it
    /// is scaled by.
    /// @param coulomb_vectors The B vectors the Coulomb is closed with, which are
    /// the plain ones whatever the exchange is made of.
    /// @note Both unrestricted entries above are written in terms of this, so the
    /// plain case is the same code as the range separated one rather than a copy of
    /// it which drifts.
    auto _unrestricted(const std::vector<std::pair<const CSparseTensor *, double>> &operators,
                       const CSparseTensor                                        &coulomb_vectors,
                       const CMolecularBasis                                      &basis,
                       const CMolecularBasis                                      &aux_basis,
                       const CPackedMatrix                                        &left_alpha,
                       const std::vector<CPackedMatrix>                           &rights_alpha,
                       const std::vector<CPackedMatrix>                           &transposed_rights_alpha,
                       const CPackedMatrix                                        &left_beta,
                       const std::vector<CPackedMatrix>                           &rights_beta,
                       const std::vector<CPackedMatrix>                           &transposed_rights_beta) const
        -> std::pair<std::vector<CPackedMatrix>, std::vector<CPackedMatrix>>;

    /// @brief Checks the factors are of one basis and of one rank.
    auto _check_factors(const CMolecularBasis            &basis,
                        const CPackedMatrix              &left,
                        const std::vector<CPackedMatrix> &rights) const -> void;

    /// @brief The memory the driver may hold when none is named.
    static constexpr size_t _default_budget = size_t{4} * 1024 * 1024 * 1024;

    /// @brief The smallest batch of auxiliary functions taken at once.
    static constexpr size_t _min_batch = 16;

    /// @brief The screening threshold of the integrals.
    double _threshold;

    /// @brief The target number of atom pairs of a block.
    size_t _block_size;

    /// @brief The memory the driver may hold, in bytes.
    size_t _budget;

    /// @brief The driver which transforms the B vectors.
    CSimdRIFockDriver _drv;
};

#endif /* SimdRIJKResponseDriver_hpp */
