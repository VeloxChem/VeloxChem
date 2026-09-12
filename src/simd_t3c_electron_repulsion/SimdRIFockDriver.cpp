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


#include "SimdRIFockDriver.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
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
CSimdRIFockDriver::compute_bq_vectors(const CMolecule        &molecule,
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

        // NOTE: every basis function on the auxiliary side of the B vectors is a
        // sum over every one of them in the integrals, so the sum over the input
        // reaches across all of the atom basis groups. Taking one group at a time
        // would make the depth of the product the functions of that group alone,
        // which is forty on average and leaves the matrix unit idle. The integrals
        // of all the groups are gathered into one buffer first, so that the depth
        // is the whole auxiliary basis and one product replaces eighty one.

        const auto nblocks = static_cast<int>(last - first);

#pragma omp parallel
        {
            std::vector<double> gathered;

#pragma omp for schedule(dynamic)
            for (int b = 0; b < nblocks; b++)
            {
                const auto iab = first + static_cast<size_t>(b);

                const auto &a_list = basis_indices[static_cast<size_t>(ab_blocks[iab].bra_index())];

                const auto &b_list = basis_indices[static_cast<size_t>(ab_blocks[iab].ket_index())];

                // NOTE: the widest the gathered buffer has to be is every atom pair
                // of the block, as no combination of it keeps more than all of
                // them. The diagonal pairs lead the off-diagonal ones and are
                // counted apart from them.

                const auto width =
                    ab_blocks[iab].number_of_diagonal_atoms() + ab_blocks[iab].number_of_pairs();

                if (width == 0) continue;

                gathered.resize(ncols * width);

                for (const auto [la, ia] : a_list)
                {
                    for (const auto [lb, jb] : b_list)
                    {
                        const auto ncomps = static_cast<size_t>((2 * la + 1) * (2 * lb + 1));

                        for (size_t m = 0; m < ncomps; m++)
                        {
                            // NOTE: a group whose block is absent, and the atom pairs
                            // a group keeps fewer of than the widest, leave zeros.
                            // The pairs are the leading ones of one ordered list, so
                            // a zero beyond the last one a group keeps adds nothing,
                            // which is what the sum of the group would have added.

                            for (const auto &in_function : in_functions)
                            {
                                auto *rows = gathered.data() + in_function.offset * width;

                                const auto iblock =
                                    batch_map[static_cast<size_t>(b) * nin_groups + in_function.group];

                                size_t npairs_in = 0;

                                const double *in_values = nullptr;

                                if (iblock != npos)
                                {
                                    npairs_in = integrals.block(iblock).number_of_pairs(
                                        la, ia, lb, jb, in_function.momentum, in_function.index);

                                    if (npairs_in > 0)
                                    {
                                        in_values = integrals.values(iblock, la, ia, lb, jb, in_function.momentum,
                                                                     in_function.index) +
                                                    m * in_function.count * npairs_in;
                                    }
                                }

                                const auto kept = std::min(npairs_in, width);

                                for (size_t r = 0; r < in_function.count; r++)
                                {
                                    auto *row = rows + r * width;

                                    if (kept > 0) std::copy(in_values + r * npairs_in, in_values + r * npairs_in + kept, row);

                                    std::fill(row + kept, row + width, 0.0);
                                }
                            }

                            for (const auto &out_function : out_functions)
                            {
                                const auto oblock = out_map[iab * nout_groups + out_function.group];

                                if (oblock == npos) continue;

                                const auto npairs_out = bq_vectors.block(oblock).number_of_pairs(
                                    la, ia, lb, jb, out_function.momentum, out_function.index);

                                if (npairs_out == 0) continue;

                                auto *out_values = bq_vectors.values(oblock, la, ia, lb, jb, out_function.momentum,
                                                                     out_function.index) +
                                                   m * out_function.count * npairs_out;

                                _matrix_product(out_function.count,
                                                std::min(npairs_out, width),
                                                ncols,
                                                1.0,
                                                metric.get() + out_function.offset * ncols,
                                                ncols,
                                                gathered.data(),
                                                width,
                                                1.0,
                                                out_values,
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
CSimdRIFockDriver::compute_y_vector(const CSparseTensor   &bq_vectors,
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
CSimdRIFockDriver::compute_fock_matrix(const CSparseTensor       &bq_vectors,
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
CSimdRIFockDriver::compute_fock_matrix(const CSparseTensor   &bq_vectors,
                                        const CMolecularBasis &basis,
                                        const CMolecularBasis &aux_basis,
                                        const CPackedMatrix   &density) const -> CPackedMatrix
{
    const auto yvector = compute_y_vector(bq_vectors, basis, aux_basis, density);

    return compute_fock_matrix(bq_vectors, basis, aux_basis, yvector);
}

namespace {  // anonymous namespace

/// @brief One auxiliary basis function of one block of the B vectors, and where
/// its values sit in that block.
struct TAuxEntry
{
    /// @brief The index of the block.
    size_t block;

    /// @brief The angular momentum of the basis function on the auxiliary side.
    int momentum;

    /// @brief The index of the basis function within its angular momentum.
    size_t index;

    /// @brief The row of the values of a combination, i.e. the angular component
    /// and the atom of the auxiliary function among those of the combination.
    size_t row;
};

}  // anonymous namespace

auto
CSimdRIFockDriver::compute_w_vectors(const CSparseTensor        &bq_vectors,
                                      const CMolecularBasis      &basis,
                                      const CMolecularBasis      &aux_basis,
                                      const CPackedMatrix        &coefficients,
                                      const size_t                qfirst,
                                      const size_t                qlast,
                                      std::vector<CPackedMatrix> &w_vectors,
                                      const bool                  accumulate) const -> void
{
    const auto nao = basis.dimensions_of_basis();

    const auto naux = aux_basis.dimensions_of_basis();

    errors::assertMsgCritical(coefficients.get_type() == mat_t::general,
                              std::string("RIJFockDriver: The orbital coefficients must be a general matrix"));

    errors::assertMsgCritical(coefficients.number_of_rows() == nao,
                              std::string("RIJFockDriver: The orbital coefficients do not match the molecular basis"));

    errors::assertMsgCritical((qfirst <= qlast) && (qlast <= naux),
                              std::string("RIJFockDriver: The range of the auxiliary basis is out of range"));

    const auto nocc = coefficients.number_of_columns();

    const auto nrange = qlast - qfirst;

    errors::assertMsgCritical(w_vectors.size() == nrange,
                              std::string("RIJFockDriver: The W matrices do not match the range of the auxiliary basis"));

    for (auto &wmat : w_vectors)
    {
        errors::assertMsgCritical((wmat.get_type() == mat_t::general) && (wmat.number_of_rows() == nao) &&
                                      (wmat.number_of_columns() == nocc),
                                  std::string("RIJFockDriver: The W matrices do not match the basis and the orbitals"));

        if (!accumulate) wmat.zero();
    }

    if ((nrange == 0) || (nocc == 0)) return;

    const auto indices = denseidx::index_functions(basis);

    const auto starts = denseidx::make_dense_starts(basis);

    const auto strides = denseidx::make_dense_strides(basis);

    const auto nmoms = static_cast<size_t>(basis.max_angular_momentum() + 1);

    const auto aux_indices = denseidx::index_functions(aux_basis);

    const auto aux_starts = denseidx::make_dense_starts(aux_basis);

    const auto aux_strides = denseidx::make_dense_strides(aux_basis);

    const auto aux_nmoms = static_cast<size_t>(aux_basis.max_angular_momentum() + 1);

    // NOTE: the blocks are gathered by the auxiliary basis function they carry, so
    // that the work is divided over the functions rather than over the blocks. A
    // block adds into the rows of every atom it touches, and two blocks of one
    // atom basis group on the auxiliary side add into the same rows, so dividing
    // over the blocks would have the threads writing over one another. Dividing
    // over the auxiliary functions gives each thread a matrix of its own.

    std::vector<std::vector<TAuxEntry>> entries(nrange);

    const auto nblocks = bq_vectors.number_of_blocks();

    for (size_t ib = 0; ib < nblocks; ib++)
    {
        const auto &block = bq_vectors.block(ib);

        const auto &c_atoms = block.c_atoms();

        const auto natoms = c_atoms.size();

        if ((natoms == 0) || (block.a_atoms().empty())) continue;

        for (const auto [lc, kc] : aux_indices[static_cast<size_t>(block.c_index())])
        {
            const auto ncomps_c = static_cast<size_t>(2 * lc + 1);

            const auto lval_c = static_cast<size_t>(lc);

            for (size_t mc = 0; mc < ncomps_c; mc++)
            {
                for (size_t n = 0; n < natoms; n++)
                {
                    const auto gq = aux_starts[static_cast<size_t>(c_atoms[n]) * aux_nmoms + lval_c] + kc +
                                    mc * aux_strides[lval_c];

                    if ((gq >= qfirst) && (gq < qlast))
                    {
                        entries[gq - qfirst].push_back(TAuxEntry{ib, lc, kc, mc * natoms + n});
                    }
                }
            }
        }
    }

    const auto *cvalues = coefficients.data();

    const auto ntasks = static_cast<int>(nrange);

    // NOTE: a threshold of zero or less takes the product always and one above one
    // takes the sum always, and neither needs the density, which costs a walk over
    // the combinations of every block to answer.

    auto use_dense = (_dense_threshold <= 0.0);

    if ((_dense_threshold > 0.0) && (_dense_threshold <= 1.0))
    {
        // NOTE: an off-diagonal pair of atoms is held once and fills two places of
        // the square, and a diagonal pair is held with the basis functions of both
        // sides and fills one. Counting the values twice over would put the density
        // above one, which a fraction cannot be.

        size_t filled = 0;

        for (size_t ib = 0; ib < bq_vectors.number_of_blocks(); ib++)
        {
            const auto &block = bq_vectors.block(ib);

            const auto &a_atoms = block.a_atoms();

            const auto &b_atoms = block.b_atoms();

            const auto natoms = block.c_atoms().size();

            if ((a_atoms.empty()) || (natoms == 0)) continue;

            size_t ndiag = 0;

            while ((ndiag < a_atoms.size()) && (a_atoms[ndiag] == b_atoms[ndiag])) ndiag++;

            for (const auto [la, ia] : indices[static_cast<size_t>(block.a_index())])
            {
                for (const auto [lb, jb] : indices[static_cast<size_t>(block.b_index())])
                {
                    for (const auto [lc, kc] : aux_indices[static_cast<size_t>(block.c_index())])
                    {
                        const auto npairs = block.number_of_pairs(la, ia, lb, jb, lc, kc);

                        if (npairs == 0) continue;

                        const auto ncomps = static_cast<size_t>((2 * la + 1) * (2 * lb + 1) * (2 * lc + 1));

                        const auto diagonal = std::min(ndiag, npairs);

                        filled += (2 * npairs - diagonal) * natoms * ncomps;
                    }
                }
            }
        }

        const auto density = (naux > 0) ? static_cast<double>(filled) / (static_cast<double>(naux) *
                                                                        static_cast<double>(nao) *
                                                                        static_cast<double>(nao))
                                        : 0.0;

        use_dense = (density >= _dense_threshold);
    }

    if (use_dense)
    {
        // NOTE: the values of one auxiliary function are scattered into a square
        // and the square is handed to a matrix product. That is more arithmetic
        // than the sum below does, and the matrix unit of the machine runs it
        // several times faster than a loop of the compiler runs the sum, which is
        // why the trade is worth taking wherever the B vectors are not sparse.

#pragma omp parallel
        {
            std::vector<double> square(nao * nao, 0.0);

#pragma omp for schedule(dynamic)
            for (int t = 0; t < ntasks; t++)
            {
                const auto iq = static_cast<size_t>(t);

                // NOTE: an auxiliary function which no block of this call carries
                // has a square of zeros, whose product adds nothing. The direct
                // mode sweeps parts of the auxiliary basis, so most of the
                // functions of a part are of that kind and skipping them is what
                // takes the half transform of a function once rather than once for
                // every part.

                if (entries[iq].empty()) continue;

                std::fill(square.begin(), square.end(), 0.0);

                for (const auto &entry : entries[iq])
                {
                    const auto &block = bq_vectors.block(entry.block);

                    const auto &a_atoms = block.a_atoms();

                    const auto &b_atoms = block.b_atoms();

                    const auto natoms = block.c_atoms().size();

                    size_t ndiag = 0;

                    while ((ndiag < a_atoms.size()) && (a_atoms[ndiag] == b_atoms[ndiag])) ndiag++;

                    const auto nq_cell = static_cast<size_t>(2 * entry.momentum + 1) * natoms;

                    for (const auto [la, ia] : indices[static_cast<size_t>(block.a_index())])
                    {
                        for (const auto [lb, jb] : indices[static_cast<size_t>(block.b_index())])
                        {
                            const auto npairs = block.number_of_pairs(la, ia, lb, jb, entry.momentum, entry.index);

                            if (npairs == 0) continue;

                            const auto ncomps_a = static_cast<size_t>(2 * la + 1);

                            const auto ncomps_b = static_cast<size_t>(2 * lb + 1);

                            const auto lval_a = static_cast<size_t>(la);

                            const auto lval_b = static_cast<size_t>(lb);

                            const auto *base =
                                bq_vectors.values(entry.block, la, ia, lb, jb, entry.momentum, entry.index);

                            for (size_t ma = 0; ma < ncomps_a; ma++)
                            {
                                for (size_t mb = 0; mb < ncomps_b; mb++)
                                {
                                    const auto *values = base + ((ma * ncomps_b + mb) * nq_cell + entry.row) * npairs;

                                    for (size_t k = 0; k < npairs; k++)
                                    {
                                        const auto irow = starts[static_cast<size_t>(a_atoms[k]) * nmoms + lval_a] +
                                                          ia + ma * strides[lval_a];

                                        const auto rrow = starts[static_cast<size_t>(b_atoms[k]) * nmoms + lval_b] +
                                                          jb + mb * strides[lval_b];

                                        square[irow * nao + rrow] = values[k];

                                        // an off-diagonal pair of atoms is held
                                        // once and fills both halves of the square

                                        if (k >= ndiag) square[rrow * nao + irow] = values[k];
                                    }
                                }
                            }
                        }
                    }
                }

                _matrix_product(nao, nocc, nao, 1.0, square.data(), nao, cvalues, nocc,
                                accumulate ? 1.0 : 0.0, w_vectors[iq].data(), nocc);
            }
        }

        return;
    }

#pragma omp parallel for schedule(dynamic) if (ntasks > 1)
    for (int t = 0; t < ntasks; t++)
    {
        const auto iq = static_cast<size_t>(t);

        auto *wvalues = w_vectors[iq].data();

        for (const auto &entry : entries[iq])
        {
            const auto &block = bq_vectors.block(entry.block);

            const auto &a_atoms = block.a_atoms();

            const auto &b_atoms = block.b_atoms();

            const auto natoms = block.c_atoms().size();

            // NOTE: the diagonal pairs of atoms lead the atom pairs of a block and
            // are carried with the basis functions of both sides, so the sum over r
            // reaches them from both. An off-diagonal pair is carried once and adds
            // into the rows of the atoms of both of its sides.

            size_t ndiag = 0;

            while ((ndiag < a_atoms.size()) && (a_atoms[ndiag] == b_atoms[ndiag])) ndiag++;

            const auto ncomps_c = static_cast<size_t>(2 * entry.momentum + 1);

            const auto nq_cell = ncomps_c * natoms;

            for (const auto [la, ia] : indices[static_cast<size_t>(block.a_index())])
            {
                for (const auto [lb, jb] : indices[static_cast<size_t>(block.b_index())])
                {
                    const auto npairs = block.number_of_pairs(la, ia, lb, jb, entry.momentum, entry.index);

                    if (npairs == 0) continue;

                    const auto ncomps_a = static_cast<size_t>(2 * la + 1);

                    const auto ncomps_b = static_cast<size_t>(2 * lb + 1);

                    const auto lval_a = static_cast<size_t>(la);

                    const auto lval_b = static_cast<size_t>(lb);

                    const auto *base = bq_vectors.values(entry.block, la, ia, lb, jb, entry.momentum, entry.index);

                    for (size_t ma = 0; ma < ncomps_a; ma++)
                    {
                        for (size_t mb = 0; mb < ncomps_b; mb++)
                        {
                            const auto *values = base + ((ma * ncomps_b + mb) * nq_cell + entry.row) * npairs;

                            for (size_t k = 0; k < npairs; k++)
                            {
                                const auto factor = values[k];

                                if (factor == 0.0) continue;

                                const auto irow = starts[static_cast<size_t>(a_atoms[k]) * nmoms + lval_a] + ia +
                                                  ma * strides[lval_a];

                                const auto rrow = starts[static_cast<size_t>(b_atoms[k]) * nmoms + lval_b] + jb +
                                                  mb * strides[lval_b];

                                auto       *wrow = wvalues + irow * nocc;

                                const auto *crow = cvalues + rrow * nocc;

#pragma omp simd
                                for (size_t s = 0; s < nocc; s++)
                                {
                                    wrow[s] += factor * crow[s];
                                }

                                if (k < ndiag) continue;

                                auto       *wrow_t = wvalues + rrow * nocc;

                                const auto *crow_t = cvalues + irow * nocc;

#pragma omp simd
                                for (size_t s = 0; s < nocc; s++)
                                {
                                    wrow_t[s] += factor * crow_t[s];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

auto
CSimdRIFockDriver::compute_w_vectors(const CSparseTensor   &bq_vectors,
                                      const CMolecularBasis &basis,
                                      const CMolecularBasis &aux_basis,
                                      const CPackedMatrix   &coefficients,
                                      const size_t           qfirst,
                                      const size_t           qlast) const -> std::vector<CPackedMatrix>
{
    errors::assertMsgCritical(qfirst <= qlast,
                              std::string("RIJFockDriver: The range of the auxiliary basis is out of range"));

    const auto nao = basis.dimensions_of_basis();

    const auto nocc = coefficients.number_of_columns();

    std::vector<CPackedMatrix> w_vectors;

    w_vectors.reserve(qlast - qfirst);

    for (size_t q = qfirst; q < qlast; q++)
    {
        w_vectors.emplace_back(nao, nocc, mat_t::general);
    }

    compute_w_vectors(bq_vectors, basis, aux_basis, coefficients, qfirst, qlast, w_vectors);

    return w_vectors;
}

auto
CSimdRIFockDriver::compute_exchange_matrix(const std::vector<CPackedMatrix> &w_vectors,
                                            CPackedMatrix                    &matrix,
                                            const double                      factor) const -> void
{
    if (w_vectors.empty()) return;

    errors::assertMsgCritical(matrix.get_type() == mat_t::symmetric,
                              std::string("RIJFockDriver: The exchange is added to a symmetric matrix"));

    const auto nao = matrix.number_of_rows();

    const auto nocc = w_vectors.front().number_of_columns();

    for (const auto &wmat : w_vectors)
    {
        errors::assertMsgCritical((wmat.get_type() == mat_t::general) && (wmat.number_of_rows() == nao) &&
                                      (wmat.number_of_columns() == nocc),
                                  std::string("RIJFockDriver: The W matrices do not match the matrix of the exchange"));
    }

    if (nocc == 0) return;

    // NOTE: the update of one auxiliary function writes the whole triangle of the
    // basis, so the functions cannot be divided over the threads without dividing
    // what they write as well. Each thread is therefore given a triangle of its
    // own and a share of the functions, and the triangles are summed into the
    // matrix at the end. The library is left to run the update on one thread: the
    // triangle is the square of the basis and holds too few blocks to spread a
    // single update over many cores, where the auxiliary functions are thousands
    // and spread perfectly.

    // NOTE: the W matrices are row major, and the column major matrix of a row
    // major array is its transpose, so the array of W is W transposed and the
    // update of it transposed times itself is W times W transposed. The upper
    // triangle of the library is the lower triangle of the array, which is the
    // triangle the packed matrix stores.

    // NOTE: the functions are still gathered into chunks before the update, so that
    // one update reads the triangle once for the whole chunk rather than once for
    // every function. The chunk is chosen from what a thread may stage, and the
    // number of triangles from what they may hold together.

    const auto row_bytes = nocc * nao * sizeof(double);

    const auto nthreads = static_cast<size_t>(omp::get_number_of_threads());

    // NOTE: the chunk is the smaller of what a thread may stage and what leaves a
    // chunk for every thread. Taking the staging alone would give a caller with few
    // W matrices, as the way which holds the B vectors has, fewer chunks than
    // threads and leave most of them with nothing to do.

    const auto staged_limit = std::max(size_t{1}, _syrk_staging / std::max(row_bytes, size_t{1}));

    const auto shared_limit = std::max(size_t{1}, (w_vectors.size() + nthreads - 1) / nthreads);

    const auto nchunk = std::max(size_t{1}, std::min({w_vectors.size(), staged_limit, shared_limit}));

    const auto nchunks = (w_vectors.size() + nchunk - 1) / nchunk;

    const auto square_bytes = nao * nao * sizeof(double);

    const auto ntriangles =
        std::max(size_t{1},
                 std::min({nthreads, nchunks, _syrk_triangles / std::max(square_bytes, size_t{1})}));

    std::vector<std::vector<double>> triangles(ntriangles);

    double staging_time = 0.0, update_time = 0.0;

    const auto profiled = (std::getenv("VLX_RIJK_PROFILE") != nullptr);

    const auto ntasks = static_cast<int>(nchunks);

#pragma omp parallel num_threads(static_cast<int>(ntriangles)) reduction(+ : staging_time, update_time)
    {
        const auto mine = static_cast<size_t>(omp_get_thread_num());

        auto &triangle = triangles[mine];

        triangle.assign(nao * nao, 0.0);

        // NOTE: every value the update reads is written by the staging below before
        // it is read, so the buffer is left with the content of its allocation
        // rather than zeroed first.

        auto staged = std::make_unique_for_overwrite<double[]>(nchunk * nocc * nao);

#pragma omp for schedule(dynamic)
        for (int itask = 0; itask < ntasks; itask++)
        {
            const auto first = static_cast<size_t>(itask) * nchunk;

            const auto count = std::min(nchunk, w_vectors.size() - first);

            const auto depth = count * nocc;

            const auto mark_staging = profiled ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point{};

            for (size_t j = 0; j < count; j++)
            {
                const auto *values = w_vectors[first + j].data();

                for (size_t irow = 0; irow < nao; irow++)
                {
                    std::copy(values + irow * nocc, values + (irow + 1) * nocc, staged.get() + irow * depth + j * nocc);
                }
            }

            if (profiled)
            {
                staging_time += std::chrono::duration<double>(std::chrono::steady_clock::now() - mark_staging).count();
            }

            const auto mark_update = profiled ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point{};

#ifdef VLX_USE_MATHLIB

            const char uplo = 'U';

            const char trans = 'T';

            auto n_arg = static_cast<lapack_int_t>(nao);

            auto k_arg = static_cast<lapack_int_t>(depth);

            auto lda = static_cast<lapack_int_t>(depth);

            auto ldc = static_cast<lapack_int_t>(nao);

            const double one = 1.0;

            dsyrk_(&uplo, &trans, &n_arg, &k_arg, &one, staged.get(), &lda, &one, triangle.data(), &ldc);

#else

            using RowMajorMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

            Eigen::Map<const RowMajorMatrix> wmap(staged.get(), static_cast<Eigen::Index>(nao),
                                                  static_cast<Eigen::Index>(depth));

            Eigen::Map<RowMajorMatrix> cmap(triangle.data(), static_cast<Eigen::Index>(nao),
                                            static_cast<Eigen::Index>(nao));

            cmap.template selfadjointView<Eigen::Lower>().rankUpdate(wmap, 1.0);

#endif /* VLX_USE_MATHLIB */

            if (profiled)
            {
                update_time += std::chrono::duration<double>(std::chrono::steady_clock::now() - mark_update).count();
            }
        }
    }

    // NOTE: only the lower triangle of each dense matrix has been written, and it
    // is the triangle the packed matrix holds, so the triangles are summed into it
    // row by row.

    auto *values = matrix.data();

    const auto nrows = static_cast<int>(nao);

#pragma omp parallel for schedule(static, 1) if (nrows > 1)
    for (int i = 0; i < nrows; i++)
    {
        const auto irow = static_cast<size_t>(i);

        auto *packed = values + irow * (irow + 1) / 2;

        for (const auto &triangle : triangles)
        {
            const auto *row = triangle.data() + irow * nao;

            for (size_t j = 0; j <= irow; j++)
            {
                packed[j] += factor * row[j];
            }
        }
    }

    if (profiled)
    {
        const auto share = static_cast<double>(ntriangles);

        std::printf("RIJK   exchange on %zu triangles, chunk %zu: staging %.3f s, update %.3f s\n", ntriangles, nchunk,
                    staging_time / share, update_time / share);

        std::fflush(stdout);
    }
}

auto
CSimdRIFockDriver::set_dense_threshold(const double threshold) -> void
{
    _dense_threshold = threshold;
}

auto
CSimdRIFockDriver::get_dense_threshold() const -> double
{
    return _dense_threshold;
}
