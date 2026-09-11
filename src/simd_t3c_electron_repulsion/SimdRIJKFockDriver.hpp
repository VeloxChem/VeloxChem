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


#ifndef SimdRIJKFockDriver_hpp
#define SimdRIJKFockDriver_hpp

#include <cstddef>
#include <vector>

#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "PackedMatrix.hpp"
#include "SimdRIFockDriver.hpp"
#include "SparseTensor.hpp"
#include "TripleSparsityPattern.hpp"

/// @brief Class CSimdRIJKFockDriver builds the Fock matrices of the resolution of
/// the identity for one molecule and one pair of bases, one matrix per call.
///
/// @note The driver holds what does not change between the calls: the inverted
/// Cholesky factor of the metric of the fitting basis and the B vectors. Both are
/// formed once by prepare. A call then contracts a density into the Y vector and
/// the Coulomb matrix, and transforms the B vectors with the orbitals into the W
/// matrices and adds their exchange.
///
/// @note The W matrices are formed for a range of the auxiliary basis at a time,
/// their exchange is added, and the same storage is reused for the next range.
/// They are rebuilt on every call, as the orbitals change, so holding all of them
/// would cost the memory of the whole auxiliary basis and save nothing. Only the
/// B vectors have to be resident, which is what the memory check is about.
///
/// @note The exchange is added with a factor the caller passes, so that a hybrid
/// functional scales it by its fraction of exact exchange. A pure functional is
/// served by a driver which never forms the B vectors and is not this one.
/// @brief How the driver forms the Fock matrices.
/// rimode::in_memory - the B vectors are formed once and held
/// rimode::direct - the integrals are formed again on every call
enum class rimode
{
    automatic,
    in_memory,
    direct
};

class CSimdRIJKFockDriver
{
   public:
    /// @brief The default constructor.
    CSimdRIJKFockDriver() = default;

    /// @brief Gets the memory the driver holds for a molecule and its bases.
    /// @param molecule The molecule to compute the Fock matrices of.
    /// @param basis The molecular basis.
    /// @param aux_basis The auxiliary molecular basis.
    /// @param threshold The screening threshold.
    /// @return The memory of the B vectors in bytes.
    /// @note The sparsity pattern of the B vectors is described to answer this,
    /// which is what the driver would do anyway and is a small part of forming
    /// them, so the answer is the memory they will take rather than an estimate of
    /// it. The W matrices are not counted, as one range of them is held at a time
    /// and is small beside the B vectors.
    auto required_memory(const CMolecule       &molecule,
                         const CMolecularBasis &basis,
                         const CMolecularBasis &aux_basis,
                         const double           threshold) const -> size_t;

    /// @brief Forms the inverted factor of the metric and the B vectors.
    /// @param molecule The molecule to compute the Fock matrices of.
    /// @param basis The molecular basis.
    /// @param aux_basis The auxiliary molecular basis.
    /// @param threshold The screening threshold.
    /// @param memory_budget The memory the driver may hold, in bytes.
    /// @note The memory is checked against the budget before the integrals are
    /// computed, from the sparsity pattern alone, so a calculation which cannot fit
    /// is told so rather than dying in the allocator after minutes of work.
    /// @param metric_threshold The eigenvalues of the metric at or below which a
    /// direction is dropped, when the metric is inverted through its square root.
    /// @param use_inverse_square_root True to invert the square root of the metric
    /// rather than its Cholesky factor.
    /// @note The Cholesky factor is the cheaper of the two by an order of
    /// magnitude and is tried first. A fitting basis which is close to linearly
    /// dependent has no Cholesky factor to invert, and the square root is inverted
    /// instead, with a warning. Setting the flag takes that way from the start.
    /// @param mode Which way the Fock matrices are formed, or automatic to hold the
    /// B vectors when they fit in the budget and to form the integrals again on
    /// every call when they do not.
    /// @note The B vectors of a large molecule do not fit in the memory of any one
    /// machine, and the direct mode is what makes such a molecule reachable. It
    /// forms the integrals once for every batch of occupied orbitals and once more
    /// for the Coulomb matrix, so it is several times slower for each Fock matrix
    /// and asks for a hundredth of the memory. The automatic choice takes the held
    /// form wherever it fits.
    auto prepare(const CMolecule       &molecule,
                 const CMolecularBasis &basis,
                 const CMolecularBasis &aux_basis,
                 const double           threshold,
                 const size_t           memory_budget,
                 const double           metric_threshold        = 1.0e-12,
                 const bool             use_inverse_square_root = false,
                 const rimode           mode                    = rimode::automatic) -> void;

    /// @brief Computes the Fock matrix of a density and a set of orbitals.
    /// @param density The density matrix, in the packed format, symmetric for a
    /// self consistent field calculation and general for a response one.
    /// @param coefficients The molecular orbital coefficients of the occupied
    /// orbitals, as a general matrix of one row per basis function and one column
    /// per orbital.
    /// @param exchange_scaling_factor The factor the exchange is scaled by, which
    /// is one for the exchange of a field calculation and the fraction of exact
    /// exchange of a hybrid functional. The exchange is not formed at all when it
    /// is zero.
    /// @return The Fock matrix, in the packed format as a symmetric matrix.
    /// @note The matrix is twice the Coulomb matrix less the scaled exchange,
    /// which is the convention of a closed shell calculation, whose density is
    /// that of one spin.
    auto compute(const CPackedMatrix &density,
                 const CPackedMatrix &coefficients,
                 const double         exchange_scaling_factor) -> CPackedMatrix;

    /// @brief Checks that the driver has been prepared.
    /// @return True if the driver is ready to form a Fock matrix.
    auto is_prepared() const -> bool;

    /// @brief Gets the way the driver forms the Fock matrices.
    /// @return The mode, which is never automatic once the driver is prepared.
    auto get_mode() const -> rimode;

    /// @brief Gets the B vectors the driver holds.
    /// @return The B vectors.
    auto get_bq_vectors() const -> const CSparseTensor &;

    /// @brief Gets the inverted Cholesky factor of the metric the driver holds.
    /// @return The inverted factor.
    auto get_metric() const -> const CPackedMatrix &;

   private:
    /// @brief Computes the Fock matrix by forming the integrals again on every
    /// call, holding no B vectors.
    /// @param density The density matrix.
    /// @param coefficients The molecular orbital coefficients.
    /// @param exchange_scaling_factor The factor the exchange is scaled by.
    /// @return The Fock matrix, twice the Coulomb less the scaled exchange.
    auto _compute_direct(const CPackedMatrix &density,
                         const CPackedMatrix &coefficients,
                         const double         exchange_scaling_factor) -> CPackedMatrix;

    /// @brief Solves the Cholesky factor of the metric against a set of right hand
    /// sides, in place.
    /// @param values The right hand sides, as a row major array of one row per
    /// auxiliary basis function and ncols columns, overwritten by the solution.
    /// @param nrows The number of auxiliary basis functions.
    /// @param ncols The number of right hand sides.
    /// @param transposed True to solve against the transpose of the factor.
    /// @note Solving rather than multiplying by an inverse, which is both cheaper
    /// and better behaved, and is why the direct way keeps the factor itself.
    auto _solve_factor(double *values, const size_t nrows, const size_t ncols, const bool transposed) const -> void;

    /// @brief The number of auxiliary basis functions whose W matrices are formed
    /// at a time.
    /// @note The exchange of a range is added before the next is formed, so this
    /// sets the memory of the W matrices rather than the work, which does not
    /// depend on it. It is the depth the rank k update of the exchange is given,
    /// which wants to be large enough to fill the cores.
    static constexpr size_t _w_batch = 64;

    /// @brief The memory a batch of the half transformed integrals is allowed to
    /// reach in the direct mode.
    /// @note The half transformed integrals are the auxiliary basis by the basis
    /// functions by the orbitals of a batch, and the batch of orbitals is chosen
    /// from this. Every batch forms the integrals again, so a larger batch is
    /// fewer passes over them and more memory.
    static constexpr size_t _direct_budget = size_t{4} * 1024 * 1024 * 1024;

    /// @brief The way the driver forms the Fock matrices.
    rimode _mode = rimode::automatic;

    /// @brief The molecule, which the direct mode forms the integrals of again on
    /// every call.
    CMolecule _molecule;

    /// @brief The sparsity pattern of the integrals, described once.
    CTripleSparsityPattern _pattern;

    /// @brief The lower triangular Cholesky factor of the metric, which the direct
    /// mode solves with.
    CPackedMatrix _factor;

    /// @brief The molecular basis.
    CMolecularBasis _basis;

    /// @brief The auxiliary molecular basis.
    CMolecularBasis _aux_basis;

    /// @brief The inverted Cholesky factor of the metric of the auxiliary basis.
    CPackedMatrix _metric;

    /// @brief The B vectors.
    CSparseTensor _bq_vectors;

    /// @brief The W matrices of one range of the auxiliary basis.
    std::vector<CPackedMatrix> _w_vectors;

    /// @brief The driver of the B vectors and of the matrices formed from them.
    CSimdRIFockDriver _drv;

    /// @brief Whether the B vectors have been formed.
    bool _prepared = false;
};

#endif /* SimdRIJKFockDriver_hpp */
