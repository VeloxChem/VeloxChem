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


#include "SimdRIJFockDriver.hpp"

#include <algorithm>
#include <limits>
#include <memory>
#include <numeric>
#include <utility>
#include <vector>

#include "AtomBasisGroup.hpp"
#include "AtomBasisTripleSparsity.hpp"
#include "DenseIndexFunc.hpp"
#include "ErrorHandler.hpp"
#include "OpenMPFunc.hpp"
#include "ScreeningFunc.hpp"
#include "SimdT3CDistributor.hpp"
#include "SimdThreeCenterElectronRepulsionDriver.hpp"
#include "TripleSparsityPattern.hpp"

#ifdef VLX_USE_MATHLIB
#include "MathLibrary.hpp"
#else
#include "Eigen/Dense"
#endif

namespace {  // anonymous namespace

/// @brief One basis function of one atom basis group on the auxiliary side, and
/// the place its auxiliary functions occupy in the permuted metric.
/// @note The values of a combination of basis functions are laid out as the
/// angular component, then the atom on the auxiliary side, then the atom pair,
/// so the auxiliary functions of a combination are the pairs of an angular
/// component and an atom, in that order. The permuted metric is ordered the same
/// way, so that the rows and the columns a contraction needs are one contiguous
/// block of it rather than a gather.
struct TAuxFunction
{
    /// @brief The index of the atom basis group on the auxiliary side.
    size_t group;

    /// @brief The angular momentum of the basis function.
    int momentum;

    /// @brief The index of the basis function within its angular momentum.
    size_t index;

    /// @brief The number of auxiliary functions, i.e. the angular components
    /// times the atoms of the group.
    size_t count;

    /// @brief The first row or column of the permuted metric.
    size_t offset;
};

/// @brief Computes the product of two row major matrices, C = alpha A B + beta C.
/// @param nrows The number of rows of A and of C.
/// @param ncols The number of columns of B and of C.
/// @param nsums The number of columns of A and of rows of B.
/// @param alpha The factor of the product.
/// @param amat The values of A, as a row major array with leading dimension lda.
/// @param lda The leading dimension of A.
/// @param bmat The values of B, as a row major array with leading dimension ldb.
/// @param ldb The leading dimension of B.
/// @param beta The factor of C.
/// @param cmat The values of C, as a row major array with leading dimension ldc.
/// @param ldc The leading dimension of C.
auto
_matrix_product(const size_t  nrows,
                const size_t  ncols,
                const size_t  nsums,
                const double  alpha,
                const double *amat,
                const size_t  lda,
                const double *bmat,
                const size_t  ldb,
                const double  beta,
                double       *cmat,
                const size_t  ldc) -> void
{
#ifdef VLX_USE_MATHLIB

    // NOTE: the library is column major, and the column major matrix of a row
    // major array is its transpose. The transpose of A B is B transposed times A
    // transposed, so the product of the row major arrays is the product of the
    // two in the other order, with the rows and the columns swapped.

    const char trans = 'N';

    auto m_arg = static_cast<lapack_int_t>(ncols);

    auto n_arg = static_cast<lapack_int_t>(nrows);

    auto k_arg = static_cast<lapack_int_t>(nsums);

    auto ldb_arg = static_cast<lapack_int_t>(ldb);

    auto lda_arg = static_cast<lapack_int_t>(lda);

    auto ldc_arg = static_cast<lapack_int_t>(ldc);

    dgemm_(&trans, &trans, &m_arg, &n_arg, &k_arg, &alpha, bmat, &ldb_arg, amat, &lda_arg, &beta, cmat, &ldc_arg);

#else

    using RowMajorMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    using RowMajorStride = Eigen::Stride<Eigen::Dynamic, 1>;

    const auto rows = static_cast<Eigen::Index>(nrows);

    const auto cols = static_cast<Eigen::Index>(ncols);

    const auto sums = static_cast<Eigen::Index>(nsums);

    Eigen::Map<const RowMajorMatrix, 0, RowMajorStride> amap(amat, rows, sums, RowMajorStride(static_cast<Eigen::Index>(lda), 1));

    Eigen::Map<const RowMajorMatrix, 0, RowMajorStride> bmap(bmat, sums, cols, RowMajorStride(static_cast<Eigen::Index>(ldb), 1));

    Eigen::Map<RowMajorMatrix, 0, RowMajorStride> cmap(cmat, rows, cols, RowMajorStride(static_cast<Eigen::Index>(ldc), 1));

    if (beta == 0.0)
    {
        cmap.noalias() = alpha * amap * bmap;
    }
    else
    {
        cmap.noalias() = beta * cmap + alpha * amap * bmap;
    }

#endif /* VLX_USE_MATHLIB */
}

/// @brief Describes the basis functions of the atom basis groups on the auxiliary
/// side and the place they occupy in the permuted metric.
/// @param groups The atom basis groups on the auxiliary side.
/// @param indices The angular momenta and indices of the basis functions of each
/// unique atom basis of the auxiliary molecular basis.
/// @param dimension The number of auxiliary functions the groups carry, which the
/// function sets.
/// @return The vector of basis functions.
auto
_make_aux_functions(const std::vector<CAtomBasisGroup>  &groups,
                    const denseidx::TBasisFunctionIndex &indices,
                    size_t                              &dimension) -> std::vector<TAuxFunction>
{
    std::vector<TAuxFunction> functions;

    size_t offset = 0;

    for (size_t i = 0; i < groups.size(); i++)
    {
        const auto natoms = groups[i].number_of_atoms();

        for (const auto [momentum, index] : indices[static_cast<size_t>(groups[i].index())])
        {
            const auto count = static_cast<size_t>(2 * momentum + 1) * natoms;

            functions.push_back(TAuxFunction{i, momentum, index, count, offset});

            offset += count;
        }
    }

    dimension = offset;

    return functions;
}

/// @brief Creates the dense indices of the auxiliary functions of the basis
/// functions, in the order the permuted metric holds them.
/// @param functions The basis functions on the auxiliary side.
/// @param groups The atom basis groups on the auxiliary side.
/// @param starts The dense index of the first angular component of the first basis
/// function of each angular momentum of each atom.
/// @param strides The distance between the dense indices of two consecutive angular
/// components of a basis function of each angular momentum.
/// @param nmoms The number of angular momenta of the auxiliary molecular basis.
/// @param dimension The number of auxiliary functions the groups carry.
/// @return The vector of dense indices.
auto
_make_dense_indices(const std::vector<TAuxFunction>    &functions,
                    const std::vector<CAtomBasisGroup> &groups,
                    const std::vector<size_t>          &starts,
                    const std::vector<size_t>          &strides,
                    const size_t                        nmoms,
                    const size_t                        dimension) -> std::vector<size_t>
{
    std::vector<size_t> dense(dimension, 0);

    for (const auto &function : functions)
    {
        const auto &atoms = groups[function.group].atoms();

        const auto natoms = atoms.size();

        const auto lval = static_cast<size_t>(function.momentum);

        const auto ncomps = static_cast<size_t>(2 * function.momentum + 1);

        for (size_t m = 0; m < ncomps; m++)
        {
            for (size_t n = 0; n < natoms; n++)
            {
                const auto atom = static_cast<size_t>(atoms[n]);

                dense[function.offset + m * natoms + n] =
                    starts[atom * nmoms + lval] + function.index + m * strides[lval];
            }
        }
    }

    return dense;
}

}  // anonymous namespace

auto
CSimdRIJFockDriver::compute_bq_vectors(const CMolecule        &molecule,
                                       const CMolecularBasis  &basis,
                                       const CMolecularBasis  &aux_basis,
                                       const CPackedMatrix    &inverse_metric,
                                       const double            threshold,
                                       const std::vector<int> &aux_atoms) const -> CSparseTensor
{
    const auto naux = aux_basis.dimensions_of_basis();

    errors::assertMsgCritical((inverse_metric.number_of_rows() == naux) && (inverse_metric.number_of_columns() == naux),
                              std::string("RIJFockDriver: The inverse metric does not match the auxiliary basis"));

    // NOTE: the metric is read as at(q, p), so the B vectors are the sum over p of
    // M_qp times the integrals of p. A symmetric matrix makes the two orders of the
    // index the same, and a lower triangular one does not: the inverted Cholesky
    // factor L, with the matrix equal to L L transposed, is what makes B transposed
    // times B the inverse of the matrix and closes the Coulomb matrix of the fitting.

    errors::assertMsgCritical((inverse_metric.get_type() == mat_t::symmetric) ||
                                  (inverse_metric.get_type() == mat_t::lower_triangular),
                              std::string("RIJFockDriver: The inverse metric must be symmetric or lower triangular"));

    // the blocks of atom pairs, which both the integrals and the B vectors carry

    auto groups = basis.basis_pair_groups();

    auto ab_blocks = sparsity::make_triple_blocks(molecule, groups, 0);

    // the atom basis groups on the auxiliary side: every one of them for the sum
    // over p, and the ones the caller asked for on q

    std::vector<int> all_atoms(static_cast<size_t>(molecule.number_of_atoms()));

    std::iota(all_atoms.begin(), all_atoms.end(), 0);

    const auto in_groups = sparsity::select_aux_groups(molecule, aux_basis, all_atoms);

    const auto out_groups = sparsity::select_aux_groups(molecule, aux_basis, aux_atoms.empty() ? all_atoms : aux_atoms);

    // the basis functions on the auxiliary side and the permuted metric

    const auto aux_indices = denseidx::index_functions(aux_basis);

    const auto aux_starts = denseidx::make_dense_starts(aux_basis);

    const auto aux_strides = denseidx::make_dense_strides(aux_basis);

    const auto nmoms = static_cast<size_t>(aux_basis.max_angular_momentum() + 1);

    size_t nrows = 0, ncols = 0;

    const auto out_functions = _make_aux_functions(out_groups, aux_indices, nrows);

    const auto in_functions = _make_aux_functions(in_groups, aux_indices, ncols);

    const auto out_dense = _make_dense_indices(out_functions, out_groups, aux_starts, aux_strides, nmoms, nrows);

    const auto in_dense = _make_dense_indices(in_functions, in_groups, aux_starts, aux_strides, nmoms, ncols);

    // NOTE: the metric is held in the order the values of a combination are laid
    // out rather than in the order of the basis, so that the rows and the columns
    // a contraction needs are a contiguous block of it. It is the square of the
    // auxiliary basis, which is nothing beside the B vectors themselves.

    auto metric = std::make_unique_for_overwrite<double[]>(nrows * ncols);

    const auto nmetric = static_cast<int>(nrows);

#pragma omp parallel for schedule(static) if (nmetric > 1)
    for (int i = 0; i < nmetric; i++)
    {
        const auto irow = static_cast<size_t>(i);

        auto *row = metric.get() + irow * ncols;

        for (size_t j = 0; j < ncols; j++)
        {
            row[j] = inverse_metric.at(out_dense[irow], in_dense[j]);
        }
    }

    // the sparsity pattern of the B vectors, and the block each pair of a block of
    // atom pairs and a basis function on the auxiliary side is held in

    const auto &bound = screenfunc::three_center_electron_repulsion_bound;

    const auto npos = std::numeric_limits<size_t>::max();

    const auto nab = ab_blocks.size();

    const auto nout_groups = out_groups.size();

    const auto nin_groups = in_groups.size();

    std::vector<CAtomBasisTripleSparsity> out_blocks;

    std::vector<size_t> out_map(nab * nout_groups, npos);

    for (size_t i = 0; i < nab; i++)
    {
        for (size_t j = 0; j < nout_groups; j++)
        {
            auto block = CAtomBasisTripleSparsity(ab_blocks[i], out_groups[j], bound, threshold);

            if (block.number_of_pairs() > 0)
            {
                out_map[i * nout_groups + j] = out_blocks.size();

                out_blocks.push_back(std::move(block));
            }
        }
    }

    auto bq_vectors = CSparseTensor(CTripleSparsityPattern(std::move(out_blocks), mat_t::symmetric, threshold));

    bq_vectors.allocate();

    bq_vectors.zero();

    // NOTE: the blocks of the integrals are described once rather than once per
    // batch, so that the batches are formed from the memory they are known to
    // need rather than from an estimate of it.

    std::vector<CAtomBasisTripleSparsity> in_blocks;

    std::vector<size_t> in_map(nab * nin_groups, npos);

    std::vector<size_t> in_elements(nab, 0);

    for (size_t i = 0; i < nab; i++)
    {
        for (size_t j = 0; j < nin_groups; j++)
        {
            auto block = CAtomBasisTripleSparsity(ab_blocks[i], in_groups[j], bound, threshold);

            if (block.number_of_pairs() > 0)
            {
                in_elements[i] += block.number_of_elements();

                in_map[i * nin_groups + j] = in_blocks.size();

                in_blocks.push_back(std::move(block));
            }
        }
    }

    const auto basis_indices = denseidx::index_functions(basis);

    CSimdThreeCenterElectronRepulsionDriver eri_drv;

    // the batches of blocks of atom pairs the integrals are formed in

    size_t first = 0;

    while (first < nab)
    {
        auto last = first;

        size_t memory = 0;

        // NOTE: a block whose integrals exceed the budget on their own is still a
        // batch of one, as the alternative is to compute nothing.

        while ((last < nab) && ((last == first) || (memory + in_elements[last] * sizeof(double) <= _batch_budget)))
        {
            memory += in_elements[last] * sizeof(double);

            last++;
        }

        std::vector<CAtomBasisTripleSparsity> batch_blocks;

        std::vector<size_t> batch_map((last - first) * nin_groups, npos);

        for (size_t i = first; i < last; i++)
        {
            for (size_t j = 0; j < nin_groups; j++)
            {
                if (const auto index = in_map[i * nin_groups + j]; index != npos)
                {
                    batch_map[(i - first) * nin_groups + j] = batch_blocks.size();

                    batch_blocks.push_back(in_blocks[index]);
                }
            }
        }

        const auto pattern = CTripleSparsityPattern(std::move(batch_blocks), mat_t::symmetric, threshold);

        auto integrals = CSparseTensor(pattern);

        integrals.allocate();

        auto distributor = CSimdT3CDistributor<CSparseTensor>(&integrals);

        eri_drv.compute(pattern, molecule, basis, aux_basis, distributor);

        // the contraction of the batch, over the blocks of atom pairs, which write
        // into blocks of the B vectors of their own and never into a shared one

        const auto nblocks = static_cast<int>(last - first);

#pragma omp parallel for schedule(dynamic) if (nblocks > 1)
        for (int b = 0; b < nblocks; b++)
        {
            const auto iab = first + static_cast<size_t>(b);

            const auto &a_list = basis_indices[static_cast<size_t>(ab_blocks[iab].bra_index())];

            const auto &b_list = basis_indices[static_cast<size_t>(ab_blocks[iab].ket_index())];

            for (const auto &out_function : out_functions)
            {
                const auto oblock = out_map[iab * nout_groups + out_function.group];

                if (oblock == npos) continue;

                const auto &out_pattern = bq_vectors.block(oblock);

                for (const auto &in_function : in_functions)
                {
                    const auto iblock = batch_map[static_cast<size_t>(b) * nin_groups + in_function.group];

                    if (iblock == npos) continue;

                    const auto &in_pattern = integrals.block(iblock);

                    const auto *submetric = metric.get() + out_function.offset * ncols + in_function.offset;

                    for (const auto [la, ia] : a_list)
                    {
                        for (const auto [lb, jb] : b_list)
                        {
                            const auto npairs_out =
                                out_pattern.number_of_pairs(la, ia, lb, jb, out_function.momentum, out_function.index);

                            const auto npairs_in =
                                in_pattern.number_of_pairs(la, ia, lb, jb, in_function.momentum, in_function.index);

                            // NOTE: the surviving atom pairs of a combination are the
                            // leading ones of the same ordered list, so the pairs the
                            // integrals carry are the pairs the B vectors carry until
                            // one of the two runs out. The bound depends on the
                            // angular momentum on the auxiliary side, so either may be
                            // the shorter, and the pairs beyond the shorter are
                            // dropped rather than written outside the block.

                            const auto npairs = std::min(npairs_in, npairs_out);

                            if (npairs == 0) continue;

                            auto *out_values =
                                bq_vectors.values(oblock, la, ia, lb, jb, out_function.momentum, out_function.index);

                            const auto *in_values =
                                integrals.values(iblock, la, ia, lb, jb, in_function.momentum, in_function.index);

                            const auto ncomps = static_cast<size_t>((2 * la + 1) * (2 * lb + 1));

                            for (size_t m = 0; m < ncomps; m++)
                            {
                                _matrix_product(out_function.count,
                                                npairs,
                                                in_function.count,
                                                1.0,
                                                submetric,
                                                ncols,
                                                in_values + m * in_function.count * npairs_in,
                                                npairs_in,
                                                1.0,
                                                out_values + m * out_function.count * npairs_out,
                                                npairs_out);
                            }
                        }
                    }
                }
            }
        }

        first = last;
    }

    return bq_vectors;
}

auto
CSimdRIJFockDriver::compute_y_vector(const CSparseTensor   &bq_vectors,
                                     const CMolecularBasis &basis,
                                     const CMolecularBasis &aux_basis,
                                     const CPackedMatrix   &density) const -> std::vector<double>
{
    const auto nao = basis.dimensions_of_basis();

    const auto naux = aux_basis.dimensions_of_basis();

    errors::assertMsgCritical((density.number_of_rows() == nao) && (density.number_of_columns() == nao),
                              std::string("RIJFockDriver: The density does not match the molecular basis"));

    errors::assertMsgCritical((density.get_type() == mat_t::symmetric) || (density.get_type() == mat_t::general),
                              std::string("RIJFockDriver: The density must be symmetric or general"));

    // NOTE: a symmetric density is stored as one triangle and its transposed
    // element is the element itself, so the off-diagonal pairs of atoms carry
    // twice it. A general density is stored in full and carries the element and
    // its transpose, which is what a response density needs.

    const auto symmetric = (density.get_type() == mat_t::symmetric);

    const auto indices = denseidx::index_functions(basis);

    const auto starts = denseidx::make_dense_starts(basis);

    const auto strides = denseidx::make_dense_strides(basis);

    const auto nmoms = static_cast<size_t>(basis.max_angular_momentum() + 1);

    const auto aux_indices = denseidx::index_functions(aux_basis);

    const auto aux_starts = denseidx::make_dense_starts(aux_basis);

    const auto aux_strides = denseidx::make_dense_strides(aux_basis);

    const auto aux_nmoms = static_cast<size_t>(aux_basis.max_angular_momentum() + 1);

    std::vector<double> yvector(naux, 0.0);

    const auto nblocks = static_cast<int>(bq_vectors.number_of_blocks());

    // NOTE: the blocks of the B vectors carry the auxiliary functions of any atom
    // of their group, so two blocks add to the same element of the Y vector. Each
    // thread sums into a vector of its own rather than through an atomic on every
    // element, and the vectors are added up in the order of the threads, so that
    // the result does not depend on the order the threads happen to finish in.

    const auto nthreads = static_cast<size_t>(omp::get_number_of_threads());

    std::vector<double> buffers(nthreads * naux, 0.0);

#pragma omp parallel
    {
        auto *partial = buffers.data() + static_cast<size_t>(omp_get_thread_num()) * naux;

        std::vector<double> weights;

        std::vector<size_t> aux_rows;

#pragma omp for schedule(dynamic)
        for (int i = 0; i < nblocks; i++)
        {
            const auto iblock = static_cast<size_t>(i);

            const auto &block = bq_vectors.block(iblock);

            const auto &a_atoms = block.a_atoms();

            const auto &b_atoms = block.b_atoms();

            const auto &c_atoms = block.c_atoms();

            const auto npairs_max = a_atoms.size();

            const auto natoms = c_atoms.size();

            if ((npairs_max == 0) || (natoms == 0)) continue;

            // NOTE: the diagonal pairs of atoms are at zero interatomic distance
            // and lead the atom pairs of a block, so counting them is enough to
            // tell the two kinds of pair apart.

            size_t ndiag = 0;

            while ((ndiag < npairs_max) && (a_atoms[ndiag] == b_atoms[ndiag])) ndiag++;

            weights.resize(npairs_max);

            const auto &a_list = indices[static_cast<size_t>(block.a_index())];

            const auto &b_list = indices[static_cast<size_t>(block.b_index())];

            const auto &c_list = aux_indices[static_cast<size_t>(block.c_index())];

            for (const auto [la, ia] : a_list)
            {
                for (const auto [lb, jb] : b_list)
                {
                    const auto ncomps_a = static_cast<size_t>(2 * la + 1);

                    const auto ncomps_b = static_cast<size_t>(2 * lb + 1);

                    const auto lval_a = static_cast<size_t>(la);

                    const auto lval_b = static_cast<size_t>(lb);

                    for (size_t ma = 0; ma < ncomps_a; ma++)
                    {
                        for (size_t mb = 0; mb < ncomps_b; mb++)
                        {
                            // NOTE: the elements of the density a combination needs
                            // do not depend on the auxiliary side, so they are
                            // gathered once and reused by every basis function of it.

                            for (size_t k = 0; k < npairs_max; k++)
                            {
                                const auto row = starts[static_cast<size_t>(a_atoms[k]) * nmoms + lval_a] + ia +
                                                 ma * strides[lval_a];

                                const auto col = starts[static_cast<size_t>(b_atoms[k]) * nmoms + lval_b] + jb +
                                                 mb * strides[lval_b];

                                if (k < ndiag)
                                {
                                    weights[k] = density.at(row, col);
                                }
                                else if (symmetric)
                                {
                                    weights[k] = 2.0 * density.at(row, col);
                                }
                                else
                                {
                                    weights[k] = density.at(row, col) + density.at(col, row);
                                }
                            }

                            for (const auto [lc, kc] : c_list)
                            {
                                const auto npairs = block.number_of_pairs(la, ia, lb, jb, lc, kc);

                                if (npairs == 0) continue;

                                const auto ncomps_c = static_cast<size_t>(2 * lc + 1);

                                const auto lval_c = static_cast<size_t>(lc);

                                const auto nq = ncomps_c * natoms;

                                const auto *values =
                                    bq_vectors.values(iblock, la, ia, lb, jb, lc, kc) + (ma * ncomps_b + mb) * nq * npairs;

                                aux_rows.resize(nq);

                                for (size_t mc = 0; mc < ncomps_c; mc++)
                                {
                                    for (size_t n = 0; n < natoms; n++)
                                    {
                                        aux_rows[mc * natoms + n] =
                                            aux_starts[static_cast<size_t>(c_atoms[n]) * aux_nmoms + lval_c] + kc +
                                            mc * aux_strides[lval_c];
                                    }
                                }

                                for (size_t q = 0; q < nq; q++)
                                {
                                    const auto *row = values + q * npairs;

                                    double sum = 0.0;

#pragma omp simd reduction(+ : sum)
                                    for (size_t p = 0; p < npairs; p++)
                                    {
                                        sum += row[p] * weights[p];
                                    }

                                    partial[aux_rows[q]] += sum;
                                }
                            }
                        }
                    }
                }
            }
        }

    }

    const auto nelems = static_cast<int>(naux);

#pragma omp parallel for schedule(static) if (nelems > 1)
    for (int k = 0; k < nelems; k++)
    {
        double sum = 0.0;

        for (size_t t = 0; t < nthreads; t++)
        {
            sum += buffers[t * naux + static_cast<size_t>(k)];
        }

        yvector[static_cast<size_t>(k)] = sum;
    }

    return yvector;
}

auto
CSimdRIJFockDriver::compute_fock_matrix(const CSparseTensor       &bq_vectors,
                                        const CMolecularBasis     &basis,
                                        const CMolecularBasis     &aux_basis,
                                        const std::vector<double> &y_vector) const -> CPackedMatrix
{
    const auto nao = basis.dimensions_of_basis();

    const auto naux = aux_basis.dimensions_of_basis();

    errors::assertMsgCritical(y_vector.size() == naux,
                              std::string("RIJFockDriver: The Y vector does not match the auxiliary basis"));

    const auto indices = denseidx::index_functions(basis);

    const auto starts = denseidx::make_dense_starts(basis);

    const auto strides = denseidx::make_dense_strides(basis);

    const auto nmoms = static_cast<size_t>(basis.max_angular_momentum() + 1);

    const auto aux_indices = denseidx::index_functions(aux_basis);

    const auto aux_starts = denseidx::make_dense_starts(aux_basis);

    const auto aux_strides = denseidx::make_dense_strides(aux_basis);

    const auto aux_nmoms = static_cast<size_t>(aux_basis.max_angular_momentum() + 1);

    auto fock = CPackedMatrix(nao, nao, mat_t::symmetric);

    fock.zero();

    auto *fock_values = fock.data();

    const auto nblocks = static_cast<int>(bq_vectors.number_of_blocks());

    // NOTE: one pair of atoms is carried by one block per atom basis group on the
    // auxiliary side, and all of them add into the same element of the Coulomb
    // matrix, so the blocks are not free of one another. Each thread sums into a
    // matrix of its own, and they are added up in the order of the threads, so
    // that the result does not depend on the order the threads finish in.

    const auto nvalues = fock.number_of_elements();

    const auto nthreads = static_cast<size_t>(omp::get_number_of_threads());

    std::vector<double> buffers(nthreads * nvalues, 0.0);

#pragma omp parallel
    {
        auto *partial = buffers.data() + static_cast<size_t>(omp_get_thread_num()) * nvalues;

        std::vector<double> contributions;

#pragma omp for schedule(dynamic)
        for (int i = 0; i < nblocks; i++)
        {
            const auto iblock = static_cast<size_t>(i);

            const auto &block = bq_vectors.block(iblock);

            const auto &a_atoms = block.a_atoms();

            const auto &b_atoms = block.b_atoms();

            const auto &c_atoms = block.c_atoms();

            const auto npairs_max = a_atoms.size();

            const auto natoms = c_atoms.size();

            if ((npairs_max == 0) || (natoms == 0)) continue;

            // NOTE: the diagonal pairs of atoms lead the atom pairs of a block, and
            // are the ones the tensor carries with the basis functions of both
            // sides, so they deliver an element of the matrix and its transpose.

            size_t ndiag = 0;

            while ((ndiag < npairs_max) && (a_atoms[ndiag] == b_atoms[ndiag])) ndiag++;

            contributions.resize(npairs_max);

            const auto &a_list = indices[static_cast<size_t>(block.a_index())];

            const auto &b_list = indices[static_cast<size_t>(block.b_index())];

            const auto &c_list = aux_indices[static_cast<size_t>(block.c_index())];

            for (const auto [la, ia] : a_list)
            {
                for (const auto [lb, jb] : b_list)
                {
                    const auto ncomps_a = static_cast<size_t>(2 * la + 1);

                    const auto ncomps_b = static_cast<size_t>(2 * lb + 1);

                    const auto lval_a = static_cast<size_t>(la);

                    const auto lval_b = static_cast<size_t>(lb);

                    for (size_t ma = 0; ma < ncomps_a; ma++)
                    {
                        for (size_t mb = 0; mb < ncomps_b; mb++)
                        {
                            std::fill(contributions.begin(), contributions.end(), 0.0);

                            auto touched = size_t{0};

                            for (const auto [lc, kc] : c_list)
                            {
                                const auto npairs = block.number_of_pairs(la, ia, lb, jb, lc, kc);

                                if (npairs == 0) continue;

                                const auto ncomps_c = static_cast<size_t>(2 * lc + 1);

                                const auto lval_c = static_cast<size_t>(lc);

                                const auto nq = ncomps_c * natoms;

                                const auto *values =
                                    bq_vectors.values(iblock, la, ia, lb, jb, lc, kc) + (ma * ncomps_b + mb) * nq * npairs;

                                touched = std::max(touched, npairs);

                                for (size_t mc = 0; mc < ncomps_c; mc++)
                                {
                                    for (size_t n = 0; n < natoms; n++)
                                    {
                                        const auto gq = aux_starts[static_cast<size_t>(c_atoms[n]) * aux_nmoms + lval_c] +
                                                        kc + mc * aux_strides[lval_c];

                                        const auto factor = y_vector[gq];

                                        const auto *row = values + (mc * natoms + n) * npairs;

#pragma omp simd
                                        for (size_t p = 0; p < npairs; p++)
                                        {
                                            contributions[p] += row[p] * factor;
                                        }
                                    }
                                }
                            }

                            // NOTE: the Coulomb matrix is symmetric and one triangle
                            // of it is stored, so an element and its transpose share
                            // a place. A diagonal pair of atoms delivers both of them
                            // and one of the two is dropped, while an off-diagonal
                            // pair delivers each unordered pair of orbitals once.

                            for (size_t k = 0; k < touched; k++)
                            {
                                const auto row = starts[static_cast<size_t>(a_atoms[k]) * nmoms + lval_a] + ia +
                                                 ma * strides[lval_a];

                                const auto col = starts[static_cast<size_t>(b_atoms[k]) * nmoms + lval_b] + jb +
                                                 mb * strides[lval_b];

                                if ((k < ndiag) && (row < col)) continue;

                                const auto upper = (row < col) ? col : row;

                                const auto lower = (row < col) ? row : col;

                                partial[upper * (upper + 1) / 2 + lower] += contributions[k];
                            }
                        }
                    }
                }
            }
        }

    }

    const auto nchunks = static_cast<int>(nvalues);

#pragma omp parallel for schedule(static) if (nchunks > 1)
    for (int k = 0; k < nchunks; k++)
    {
        double sum = 0.0;

        for (size_t t = 0; t < nthreads; t++)
        {
            sum += buffers[t * nvalues + static_cast<size_t>(k)];
        }

        fock_values[static_cast<size_t>(k)] = sum;
    }

    return fock;
}

auto
CSimdRIJFockDriver::compute_fock_matrix(const CSparseTensor   &bq_vectors,
                                        const CMolecularBasis &basis,
                                        const CMolecularBasis &aux_basis,
                                        const CPackedMatrix   &density) const -> CPackedMatrix
{
    const auto yvector = compute_y_vector(bq_vectors, basis, aux_basis, density);

    return compute_fock_matrix(bq_vectors, basis, aux_basis, yvector);
}
