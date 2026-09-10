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


#ifndef SimdRIFockDriver_hpp
#define SimdRIFockDriver_hpp

#include <cstddef>
#include <vector>

#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "PackedMatrix.hpp"
#include "SparseTensor.hpp"

/// @brief Class CSimdRIFockDriver computes the B vectors of the resolution of the
/// identity approximation of the Coulomb interaction.
///
/// @note The B vectors are defined as B(q)_ij = sum over p of (ij|p) L(-1)_pq, where
/// (ij|p) are the three-center electron repulsion integrals and L(-1) is the inverse
/// of the two-center metric (p|q), which the caller supplies in the packed format.
///
/// @note The B vectors are screened on the atomic orbital pair, as the three-center
/// integrals are, and are dense in q, as the inverse metric is dense. Their memory is
/// therefore the surviving atomic orbital pairs times the auxiliary basis, which
/// reaches hundreds of gigabytes for a large molecule. The auxiliary functions of q
/// are what a calculation distributes over the MPI ranks, and the driver forms the
/// ones of the given atoms alone, so that a rank holds its own share and no more.
/// A single rank asks for all of them.
///
/// @note The sum over p runs over the whole auxiliary basis whatever share of q a
/// rank holds, so the three-center integrals cannot be distributed the same way. They
/// are formed in batches of the blocks of atomic orbital pairs instead, and each
/// batch is contracted and dropped before the next is formed, so that the peak is the
/// B vectors and one batch rather than the whole three-center tensor.
class CSimdRIFockDriver
{
   public:
    /// @brief The default constructor.
    CSimdRIFockDriver() = default;

    /// @brief Computes the B vectors of the resolution of the identity approximation.
    /// @param molecule The molecule to compute the B vectors of.
    /// @param basis The molecular basis on a and b sides.
    /// @param aux_basis The auxiliary molecular basis.
    /// @param inverse_metric The inverse of the two-center metric, in the packed
    /// format, as packlin::invert returns it.
    /// @param threshold The screening threshold.
    /// @param aux_atoms The atoms whose auxiliary functions the B vectors are formed
    /// for, as their indices in the molecule and without repetition, or empty for all
    /// of them.
    /// @return The sparse tensor of the B vectors, whose blocks hold the auxiliary
    /// functions of q in the place the three-center integrals hold those of p.
    auto compute_bq_vectors(const CMolecule        &molecule,
                            const CMolecularBasis  &basis,
                            const CMolecularBasis  &aux_basis,
                            const CPackedMatrix    &inverse_metric,
                            const double            threshold,
                            const std::vector<int> &aux_atoms = {}) const -> CSparseTensor;

    /// @brief Contracts the B vectors with a density matrix.
    /// @param bq_vectors The B vectors, as compute_bq_vectors returns them.
    /// @param basis The molecular basis on a and b sides.
    /// @param aux_basis The auxiliary molecular basis.
    /// @param density The density matrix, in the packed format.
    /// @return The Y vector, Y(q) = sum over i and j of B(q)_ij D_ij, indexed by the
    /// dense index of the auxiliary basis function and zero where the B vectors carry
    /// no auxiliary function of their own.
    /// @note A symmetric density is the one of a self consistent field calculation
    /// and a general density is the one of a response calculation. Which of the two
    /// is contracted follows from the type of the matrix, so that the caller cannot
    /// ask for the wrong one.
    /// @note The B vectors hold each unordered pair of atoms once, and hold the
    /// diagonal pairs of atoms with the basis functions of both sides. The elements
    /// of a diagonal pair of atoms are therefore summed as they are, and those of an
    /// off-diagonal pair are summed with the transposed element of the density, which
    /// is twice the element itself when the density is symmetric.
    auto compute_y_vector(const CSparseTensor    &bq_vectors,
                          const CMolecularBasis  &basis,
                          const CMolecularBasis  &aux_basis,
                          const CPackedMatrix    &density) const -> std::vector<double>;

    /// @brief Computes the Coulomb matrix of the resolution of the identity.
    /// @param bq_vectors The B vectors, as compute_bq_vectors returns them.
    /// @param basis The molecular basis on a and b sides.
    /// @param aux_basis The auxiliary molecular basis.
    /// @param y_vector The Y vector, as compute_y_vector returns it.
    /// @return The Coulomb matrix, F_ij = sum over q of B(q)_ij Y(q), in the packed
    /// format as a symmetric matrix.
    /// @note This is the form which takes the Y vector, so that a calculation which
    /// contracts many densities with one set of B vectors forms them once.
    auto compute_fock_matrix(const CSparseTensor        &bq_vectors,
                             const CMolecularBasis      &basis,
                             const CMolecularBasis      &aux_basis,
                             const std::vector<double>  &y_vector) const -> CPackedMatrix;

    /// @brief Computes the Coulomb matrix of the resolution of the identity for a
    /// density matrix.
    /// @param bq_vectors The B vectors, as compute_bq_vectors returns them.
    /// @param basis The molecular basis on a and b sides.
    /// @param aux_basis The auxiliary molecular basis.
    /// @param density The density matrix, in the packed format.
    /// @return The Coulomb matrix, in the packed format as a symmetric matrix.
    /// @note The Coulomb matrix is symmetric whether the density is symmetric or
    /// general, so it is stored as a symmetric matrix in either case.
    /// @note No factor of the occupancy or of the spin is applied. This is the
    /// contraction as it is written, and a convention which carries such a factor
    /// applies it to the density or to the matrix returned.
    auto compute_fock_matrix(const CSparseTensor   &bq_vectors,
                             const CMolecularBasis &basis,
                             const CMolecularBasis &aux_basis,
                             const CPackedMatrix   &density) const -> CPackedMatrix;

    /// @brief Transforms one index of the B vectors into the molecular orbitals.
    /// @param bq_vectors The B vectors, as compute_bq_vectors returns them.
    /// @param basis The molecular basis on a and b sides.
    /// @param aux_basis The auxiliary molecular basis.
    /// @param coefficients The molecular orbital coefficients of the orbitals to
    /// transform into, as a general matrix of one row per basis function and one
    /// column per orbital.
    /// @param qfirst The first auxiliary basis function to transform, as its dense
    /// index in the auxiliary basis.
    /// @param qlast The auxiliary basis function past the last one to transform.
    /// @return The W matrices, W(q)_is = sum over r of B(q)_ir C_rs, one general
    /// matrix of one row per basis function and one column per orbital for each
    /// auxiliary basis function of the range.
    /// @note The sum over r leaves nothing of the sparsity of the pairs of basis
    /// functions: an atom on the auxiliary side is reached by every atom of the
    /// molecule, as the pair of an atom with itself is at zero distance and
    /// survives any threshold. The W matrices are therefore dense, and are the
    /// dimensions of the basis times the orbitals times the auxiliary basis, which
    /// grows as the cube of the size of the molecule with nothing to screen. They
    /// are formed for a range of the auxiliary basis rather than for all of it, so
    /// that the caller forms what it can hold, uses it and asks for the next range.
    auto compute_w_vectors(const CSparseTensor   &bq_vectors,
                           const CMolecularBasis &basis,
                           const CMolecularBasis &aux_basis,
                           const CPackedMatrix   &coefficients,
                           const size_t           qfirst,
                           const size_t           qlast) const -> std::vector<CPackedMatrix>;

    /// @brief Transforms one index of the B vectors into the molecular orbitals,
    /// into matrices the caller holds.
    /// @param bq_vectors The B vectors, as compute_bq_vectors returns them.
    /// @param basis The molecular basis on a and b sides.
    /// @param aux_basis The auxiliary molecular basis.
    /// @param coefficients The molecular orbital coefficients.
    /// @param qfirst The first auxiliary basis function to transform.
    /// @param qlast The auxiliary basis function past the last one to transform.
    /// @param w_vectors The matrices to store the result in, one per auxiliary
    /// basis function of the range, each of one row per basis function and one
    /// column per orbital. They are set rather than added to.
    /// @note This is the form which does not allocate, so that a calculation which
    /// forms the W matrices of one range after another reuses one set of them
    /// rather than allocating and freeing gigabytes on every range.
    auto compute_w_vectors(const CSparseTensor         &bq_vectors,
                           const CMolecularBasis       &basis,
                           const CMolecularBasis       &aux_basis,
                           const CPackedMatrix         &coefficients,
                           const size_t                 qfirst,
                           const size_t                 qlast,
                           std::vector<CPackedMatrix>  &w_vectors) const -> void;

    /// @brief Adds the exchange contribution of a range of the auxiliary basis to
    /// a matrix.
    /// @param w_vectors The W matrices of the range, as compute_w_vectors returns
    /// them.
    /// @param matrix The symmetric matrix to add the contribution to, in the packed
    /// format, which the Coulomb matrix of the same calculation has been formed in.
    /// @param factor The factor the contribution is added with.
    /// @note The contribution is factor times the sum over q of W(q) times W(q)
    /// transposed, summed over the range the W matrices hold. The ranges of a
    /// calculation are added one after another into the same matrix, as the W
    /// matrices of the whole auxiliary basis do not fit in the memory of a large
    /// molecule.
    /// @note No sign is applied. The exchange enters a Fock matrix with a sign and
    /// a factor which depend on the convention of the caller, which passes them as
    /// the factor rather than having them built in here.
    auto compute_exchange_matrix(const std::vector<CPackedMatrix> &w_vectors,
                                 CPackedMatrix                    &matrix,
                                 const double                      factor = 1.0) const -> void;

   private:
    /// @brief The number of auxiliary basis functions whose W matrices are staged
    /// into one buffer before the rank k update.
    /// @note The update of one auxiliary function alone is too small a piece of
    /// work for the library to spread over the cores, as its depth is the orbitals
    /// alone. Staging several of them makes one update of that many times the
    /// depth. Measured on tagrisso with def2-svp, in gigaflops per second: 274 for
    /// one, 442 for sixteen and 526 for sixty four, which is what the library
    /// reaches on this machine, so there is nothing further to gain. The staging
    /// buffer is this many times the orbitals times the dimensions of the basis,
    /// which is a hundred megabytes at the sizes the W matrices are formed in.
    static constexpr size_t _syrk_chunk = 64;

    /// @brief The memory a batch of the three-center integrals is allowed to reach.
    /// @note The batch is the blocks of atomic orbital pairs whose integrals are
    /// formed at once. Making it larger costs memory and buys nothing beyond the
    /// point where the threads are busy, as the contraction of a batch is what is
    /// parallelized and not the batches.
    static constexpr size_t _batch_budget = size_t{4} * 1024 * 1024 * 1024;
};

#endif /* SimdRIFockDriver_hpp */
