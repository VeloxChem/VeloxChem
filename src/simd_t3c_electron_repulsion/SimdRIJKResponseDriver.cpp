//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//

#include "SimdRIJKResponseDriver.hpp"

#include <algorithm>
#include <string>

#include "ErrorHandler.hpp"

#ifdef VLX_USE_MATHLIB
#include "MathLibrary.hpp"
#else
#include "Eigen/Dense"
#endif

/// @brief The product of a row major matrix by the transpose of another, added to
/// what is there: C += A B^T, for A of nrows by nsums and B of ncols by nsums.
/// @note The rows of the three are given as they are laid out and not as they are
/// used, so a block of a wider array is multiplied where it stands. The right
/// factors of a batch are transformed as one wide matrix, and each density's share
/// of the result is such a block.
static auto
_add_multiply_by_transpose(const size_t  nrows,
                           const size_t  ncols,
                           const size_t  nsums,
                           const double *amat,
                           const size_t  arow,
                           const double *bmat,
                           const size_t  brow,
                           double       *cmat,
                           const size_t  crow) -> void
{
#ifdef VLX_USE_MATHLIB

    const char trans_t = 'T';

    const char trans_n = 'N';

    const double alpha = 1.0;

    const double beta = 1.0;

    auto m_arg = static_cast<lapack_int_t>(ncols);

    auto n_arg = static_cast<lapack_int_t>(nrows);

    auto k_arg = static_cast<lapack_int_t>(nsums);

    auto lda_arg = static_cast<lapack_int_t>(brow);

    auto ldb_arg = static_cast<lapack_int_t>(arow);

    auto ldc_arg = static_cast<lapack_int_t>(crow);

    dgemm_(&trans_t, &trans_n, &m_arg, &n_arg, &k_arg, &alpha, bmat, &lda_arg, amat, &ldb_arg, &beta, cmat, &ldc_arg);

#else

    using RowMajorMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    using Stride = Eigen::Stride<Eigen::Dynamic, 1>;

    Eigen::Map<const RowMajorMatrix, 0, Stride> a(
        amat, static_cast<Eigen::Index>(nrows), static_cast<Eigen::Index>(nsums), Stride(static_cast<Eigen::Index>(arow), 1));

    Eigen::Map<const RowMajorMatrix, 0, Stride> b(
        bmat, static_cast<Eigen::Index>(ncols), static_cast<Eigen::Index>(nsums), Stride(static_cast<Eigen::Index>(brow), 1));

    Eigen::Map<RowMajorMatrix, 0, Stride> c(
        cmat, static_cast<Eigen::Index>(nrows), static_cast<Eigen::Index>(ncols), Stride(static_cast<Eigen::Index>(crow), 1));

    c.noalias() += a * b.transpose();

#endif
}

auto
CSimdRIJKResponseDriver::get_threshold() const -> double
{
    return _threshold;
}

auto
CSimdRIJKResponseDriver::get_block_size() const -> size_t
{
    return _block_size;
}

auto
CSimdRIJKResponseDriver::get_memory_budget() const -> size_t
{
    return _budget;
}

auto
CSimdRIJKResponseDriver::_check_factors(const CMolecularBasis            &basis,
                                        const CPackedMatrix              &left,
                                        const std::vector<CPackedMatrix> &rights) const -> void
{
    const auto nao = basis.dimensions_of_basis();

    errors::assertMsgCritical(
        left.number_of_rows() == nao,
        std::string("SimdRIJKResponseDriver: The left factor is not of the molecular basis"));

    errors::assertMsgCritical(
        left.number_of_columns() > 0,
        std::string("SimdRIJKResponseDriver: The left factor has no columns"));

    for (const auto &right : rights)
    {
        errors::assertMsgCritical(
            right.number_of_rows() == nao,
            std::string("SimdRIJKResponseDriver: A right factor is not of the molecular basis"));

        // NOTE: the two factors close over the same index, so a density whose
        // factors disagree on it is not a density. Caught here rather than by the
        // multiply, which would read past one of them.
        errors::assertMsgCritical(
            right.number_of_columns() == left.number_of_columns(),
            std::string("SimdRIJKResponseDriver: A right factor is not of the rank of the left factor"));
    }
}

auto
CSimdRIJKResponseDriver::compute_exchange(const CSparseTensor              &bq_vectors,
                                          const CMolecularBasis            &basis,
                                          const CMolecularBasis            &aux_basis,
                                          const CPackedMatrix              &left,
                                          const std::vector<CPackedMatrix> &rights) const -> std::vector<CPackedMatrix>
{
    auto exchanges = std::vector<CPackedMatrix>();

    if (rights.empty()) return exchanges;

    _check_factors(basis, left, rights);

    const auto nao = basis.dimensions_of_basis();

    const auto nvec = left.number_of_columns();

    const auto ndens = rights.size();

    const auto naux = aux_basis.dimensions_of_basis();

    for (size_t idens = 0; idens < ndens; idens++)
    {
        exchanges.push_back(CPackedMatrix(nao, nao, mat_t::general));

        exchanges.back().zero();
    }

    // NOTE: the right factors side by side as one matrix of the basis by the rank
    // times the densities. The transformation reads the B vectors, which is the
    // dearest thing it does, and one wide matrix reads them once where a matrix
    // for each density reads them once for each density.

    const auto nwide = nvec * ndens;

    auto stacked = CPackedMatrix(nao, nwide, mat_t::general);

    stacked.zero();

    for (size_t idens = 0; idens < ndens; idens++)
    {
        const auto *values = rights[idens].data();

        auto *target = stacked.data();

        for (size_t mu = 0; mu < nao; mu++)
        {
            for (size_t k = 0; k < nvec; k++)
            {
                target[mu * nwide + idens * nvec + k] = values[mu * nvec + k];
            }
        }
    }

    // NOTE: the auxiliary basis in batches, so the transformed vectors are bounded
    // by the budget and not by the problem. One batch holds the left factor
    // transformed and the right factors transformed, which is the basis by the
    // rank and the basis by the rank times the densities.

    const auto per_function = nao * nvec * (1 + ndens) * sizeof(double);

    const auto by_memory = _budget / std::max(per_function, size_t{1});

    const auto nbatch = std::min(naux, std::max(_min_batch, by_memory));

    for (size_t first = 0; first < naux; first += nbatch)
    {
        const auto last = std::min(first + nbatch, naux);

        const auto uvecs = _drv.compute_w_vectors(bq_vectors, basis, aux_basis, left, first, last);

        const auto pvecs = _drv.compute_w_vectors(bq_vectors, basis, aux_basis, stacked, first, last);

        const auto count = static_cast<int>(last - first);

        // NOTE: the densities and not the auxiliary functions are what the threads
        // divide, so that each of them writes into an exchange matrix of its own
        // and no two of them accumulate into the same one. A batch of trial vectors
        // is what a response calculation always has.

        const auto nrange = static_cast<int>(ndens);

#pragma omp parallel for schedule(static) if (nrange > 1)
        for (int at = 0; at < nrange; at++)
        {
            const auto idens = static_cast<size_t>(at);

            for (int q = 0; q < count; q++)
            {
                const auto iq = static_cast<size_t>(q);

                _add_multiply_by_transpose(nao,
                                           nao,
                                           nvec,
                                           uvecs[iq].data(),
                                           nvec,
                                           pvecs[iq].data() + idens * nvec,
                                           nwide,
                                           exchanges[idens].data(),
                                           nao);
            }
        }
    }

    return exchanges;
}
