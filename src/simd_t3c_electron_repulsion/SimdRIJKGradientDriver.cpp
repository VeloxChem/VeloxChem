//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//

#include "SimdRIJKGradientDriver.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <set>
#include <string>

#include "DenseIndexFunc.hpp"
#include "ErrorHandler.hpp"
#include "SimdRIFockCommon.hpp"
#include "SimdThreeCenterElectronRepulsionDriver.hpp"
#include "SimdThreeCenterElectronRepulsionGradientDriver.hpp"
#include "SimdThreeCenterElectronRepulsionGradientRsDriver.hpp"
#include "SimdTwoCenterElectronRepulsionGradientDriver.hpp"
#include "SimdTwoCenterElectronRepulsionGradientRsDriver.hpp"
#include "TensorComponents.hpp"

#ifdef VLX_USE_MATHLIB
#include "MathLibrary.hpp"
#else
#include "Eigen/Dense"
#endif

auto
CSimdRIJKGradientDriver::get_threshold() const -> double
{
    return _threshold;
}

auto
CSimdRIJKGradientDriver::get_block_size() const -> size_t
{
    return _block_size;
}

/// @brief The product of two row major matrices, C = A B.
/// @note The library is column major and the column major matrix of a row major
/// array is its transpose, so the product of the row major arrays is the product
/// of the two in the other order with the rows and the columns swapped. The same
/// reading the Fock driver's multiply takes.
static auto
_multiply(const size_t  nrows,
          const size_t  ncols,
          const size_t  nsums,
          const double *amat,
          const double *bmat,
          double       *cmat) -> void
{
#ifdef VLX_USE_MATHLIB

    const char trans = 'N';

    const double alpha = 1.0;

    const double beta = 0.0;

    auto m_arg = static_cast<lapack_int_t>(ncols);

    auto n_arg = static_cast<lapack_int_t>(nrows);

    auto k_arg = static_cast<lapack_int_t>(nsums);

    auto ldb_arg = static_cast<lapack_int_t>(ncols);

    auto lda_arg = static_cast<lapack_int_t>(nsums);

    auto ldc_arg = static_cast<lapack_int_t>(ncols);

    dgemm_(&trans, &trans, &m_arg, &n_arg, &k_arg, &alpha, bmat, &ldb_arg, amat, &lda_arg, &beta, cmat, &ldc_arg);

#else

    using RowMajorMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    Eigen::Map<const RowMajorMatrix> a(amat, static_cast<Eigen::Index>(nrows), static_cast<Eigen::Index>(nsums));

    Eigen::Map<const RowMajorMatrix> b(bmat, static_cast<Eigen::Index>(nsums), static_cast<Eigen::Index>(ncols));

    Eigen::Map<RowMajorMatrix> c(cmat, static_cast<Eigen::Index>(nrows), static_cast<Eigen::Index>(ncols));

    c.noalias() = a * b;

#endif
}

/// @brief The product of two row major matrices added to what is there,
/// C += A B.
/// @note The same reading of the library's ordering as the product above, with the
/// one difference that it accumulates. An open shell adds a term per spin into one
/// matrix, and the overwriting form silently kept only the last of them.
static auto
_add_multiply(const size_t  nrows,
              const size_t  ncols,
              const size_t  nsums,
              const double *amat,
              const double *bmat,
              double       *cmat) -> void
{
#ifdef VLX_USE_MATHLIB

    const char trans = 'N';

    const double alpha = 1.0;

    const double beta = 1.0;

    auto m_arg = static_cast<lapack_int_t>(ncols);

    auto n_arg = static_cast<lapack_int_t>(nrows);

    auto k_arg = static_cast<lapack_int_t>(nsums);

    auto ldb_arg = static_cast<lapack_int_t>(ncols);

    auto lda_arg = static_cast<lapack_int_t>(nsums);

    auto ldc_arg = static_cast<lapack_int_t>(ncols);

    dgemm_(&trans, &trans, &m_arg, &n_arg, &k_arg, &alpha, bmat, &ldb_arg, amat, &lda_arg, &beta, cmat, &ldc_arg);

#else

    using RowMajorMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    Eigen::Map<const RowMajorMatrix> a(amat, static_cast<Eigen::Index>(nrows), static_cast<Eigen::Index>(nsums));

    Eigen::Map<const RowMajorMatrix> b(bmat, static_cast<Eigen::Index>(nsums), static_cast<Eigen::Index>(ncols));

    Eigen::Map<RowMajorMatrix> c(cmat, static_cast<Eigen::Index>(nrows), static_cast<Eigen::Index>(ncols));

    c.noalias() += a * b;

#endif
}

/// @brief The product of a row major matrix by the transpose of another, with the
/// result added to what is there: C += A B^T, for A of nrows by nsums and B of
/// ncols by nsums.
/// @note The same reading of the library's ordering as the product above. The two
/// arrays may be the same one, which is how a Gram product is taken without a
/// second copy of the matrix transposed.
static auto
_add_multiply_by_transpose(const size_t  nrows,
                           const size_t  ncols,
                           const size_t  nsums,
                           const double *amat,
                           const double *bmat,
                           double       *cmat) -> void
{
#ifdef VLX_USE_MATHLIB

    const char trans_t = 'T';

    const char trans_n = 'N';

    const double alpha = 1.0;

    const double beta = 1.0;

    auto m_arg = static_cast<lapack_int_t>(ncols);

    auto n_arg = static_cast<lapack_int_t>(nrows);

    auto k_arg = static_cast<lapack_int_t>(nsums);

    auto lda_arg = static_cast<lapack_int_t>(nsums);

    auto ldb_arg = static_cast<lapack_int_t>(nsums);

    auto ldc_arg = static_cast<lapack_int_t>(ncols);

    dgemm_(&trans_t, &trans_n, &m_arg, &n_arg, &k_arg, &alpha, bmat, &lda_arg, amat, &ldb_arg, &beta, cmat, &ldc_arg);

#else

    using RowMajorMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    Eigen::Map<const RowMajorMatrix> a(amat, static_cast<Eigen::Index>(nrows), static_cast<Eigen::Index>(nsums));

    Eigen::Map<const RowMajorMatrix> b(bmat, static_cast<Eigen::Index>(ncols), static_cast<Eigen::Index>(nsums));

    Eigen::Map<RowMajorMatrix> c(cmat, static_cast<Eigen::Index>(nrows), static_cast<Eigen::Index>(ncols));

    c.noalias() += a * b.transpose();

#endif
}

/// @brief Applies the transpose of the inverted metric to a set of columns held
/// one after another.
/// @param triangular True where the factor is the inverted Cholesky factor, which
/// is lower triangular, and false where it is the inverted square root, which is
/// symmetric.
/// @note Both close the same sum -- M transposed times M is the inverse of the
/// metric for either of them, which is why the energy takes either. What differs is
/// the storage: a triangular factor lets the element p of the result read the rows
/// at and below p alone, and a symmetric root does not. Reading a root that way
/// drops every term below the diagonal and answers something which is not a
/// gradient of anything.
/// @note A range separated gradient cannot choose. The attenuated metric is
/// numerically singular and its Cholesky factorization fails on some fitting sets
/// and succeeds on others, so both kinds arrive and both have to work.
/// @brief The same, over a factor which the caller has already expanded.
///
/// @note The expansion is the square of the auxiliary basis and is the same array
/// every time. Forming it inside this routine costs nothing when the routine is
/// called once, and everything when it is called once per element of a matrix:
/// ninety-two megabytes written eight thousand times over on a molecule of this
/// size, which is memory traffic of a wholly different order to the arithmetic it
/// serves. The caller which loops expands it once and hands it in.
static auto
_multiply_transposed(const double *factor,
                     const size_t  naux,
                     double       *values,
                     const size_t  ncols,
                     const bool    triangular) -> void
{
    auto column = std::vector<double>(naux, 0.0);

    for (size_t icol = 0; icol < ncols; icol++)
    {
        auto *target = values + icol * naux;

        std::fill(column.begin(), column.end(), 0.0);

        for (size_t p = 0; p < naux; p++)
        {
            double sum = 0.0;

            // NOTE: a lower triangular factor has nothing above its diagonal, so
            // the element p of the result reads the rows at and below p alone. A
            // symmetric root has to be read whole.

            for (size_t q = (triangular ? p : 0); q < naux; q++)
            {
                sum += factor[q * naux + p] * target[q];
            }

            column[p] = sum;
        }

        std::copy(column.begin(), column.end(), target);
    }
}

static auto
_multiply_transposed(const CPackedMatrix &metric, double *values, const size_t ncols) -> void
{
    const auto naux = metric.number_of_rows();

    auto factor = std::vector<double>(naux * naux, 0.0);

    metric.to_dense(factor.data());

    _multiply_transposed(factor.data(), naux, values, ncols, metric.get_type() != mat_t::symmetric);
}

auto
CSimdRIJKGradientDriver::_apply_transposed_factor(const CPackedMatrix &metric,
                                                  double              *values,
                                                  const size_t         ncols) const -> void
{
    _multiply_transposed(metric, values, ncols);
}

auto
CSimdRIJKGradientDriver::_apply_transposed_factor(const CPackedMatrix        &metric,
                                                  std::vector<CPackedMatrix> &matrices) const -> void
{
    // NOTE: the auxiliary index is the one the factor is applied over, and it is
    // the index of the vector rather than an index inside a matrix. The elements
    // of one place of the matrices are therefore gathered into a column, the
    // factor applied to it, and the column scattered back.

    if (matrices.empty()) return;

    const auto naux = matrices.size();

    const auto nelements = matrices.front().number_of_elements();

    // NOTE: the transpose of the factor, written out once. For the inverted
    // Cholesky factor the element of the result at p reads the rows at and below p,
    // which is the transpose of a lower triangular matrix read as it stands. For the
    // inverted square root there is nothing to skip: it is symmetric, so its
    // transpose is itself and the whole of it is written. The expansion it is built
    // from is released before the arrays which are the size of the whole phase are
    // taken.

    const auto triangular = (metric.get_type() != mat_t::symmetric);

    auto transposed = std::vector<double>(naux * naux, 0.0);

    {
        auto factor = std::vector<double>(naux * naux, 0.0);

        metric.to_dense(factor.data());

        for (size_t p = 0; p < naux; p++)
        {
            for (size_t q = (triangular ? p : 0); q < naux; q++)
            {
                transposed[p * naux + q] = factor[q * naux + p];
            }
        }
    }

    // NOTE: the auxiliary index is the one the factor is applied over, and it is
    // the index of the vector rather than an index inside a matrix. Gathering the
    // matrices into one array indexed by the auxiliary function and then by the
    // element makes that index the one a single product contracts, so the whole
    // phase is one multiply of the library's instead of one of ours for every
    // element of a matrix -- eight thousand calls on a molecule of any size, each
    // of them a multiply of a matrix by a single column.

    // NOTE: the elements in panels, so the two arrays which are the auxiliary
    // basis times the elements are bounded by the budget rather than by the
    // problem. The factor is the square of the auxiliary basis and cannot be
    // divided, so it is taken off the budget before the panel is sized. At the
    // default budget a molecule of any size this driver holds B vectors for is
    // one panel, and the division costs nothing; it is there for the case the
    // note above describes, where the elements alone are gigabytes.

    const auto fixed = naux * naux * sizeof(double);

    const auto spare = (_budget > fixed) ? _budget - fixed : size_t{0};

    const auto per_element = 2 * naux * sizeof(double);

    const auto by_memory = spare / std::max(per_element, size_t{1});

    const auto npanel = std::min(nelements, std::max(_min_batch, by_memory));

    auto gathered = std::vector<double>(naux * npanel, 0.0);

    auto result = std::vector<double>(naux * npanel, 0.0);

    const auto nrange = static_cast<int>(naux);

    for (size_t first = 0; first < nelements; first += npanel)
    {
        const auto count = std::min(npanel, nelements - first);

        // NOTE: the panel is packed with the width it actually has and not with
        // the width of the storage, so the last and shorter one is a matrix the
        // multiply reads as it stands rather than one with a stride of its own.

#pragma omp parallel for schedule(static)
        for (int at = 0; at < nrange; at++)
        {
            const auto q = static_cast<size_t>(at);

            std::copy_n(matrices[q].data() + first, count, gathered.data() + q * count);
        }

        _multiply(naux, count, naux, transposed.data(), gathered.data(), result.data());

#pragma omp parallel for schedule(static)
        for (int at = 0; at < nrange; at++)
        {
            const auto q = static_cast<size_t>(at);

            std::copy_n(result.data() + q * count, count, matrices[q].data() + first);
        }
    }
}

auto
CSimdRIJKGradientDriver::_close_orbitals(const std::vector<double> &transposed,
                                         const size_t               nao,
                                         const size_t               norbs,
                                         const CPackedMatrix       &half) const -> CPackedMatrix
{
    // NOTE: the second half of the transformation, d(q)_ij = sum over mu of
    // C_mu,i W(q)_mu,j, which is one product of the transposed coefficients with
    // the half transformed matrix. The result is symmetric and only its lower
    // triangle is kept, but the product forms the square: the symmetric product
    // of the library would need the two factors to be the same matrix, and these
    // are not.

    auto dense_w = std::vector<double>(nao * norbs, 0.0);

    half.to_dense(dense_w.data());

    auto square = std::vector<double>(norbs * norbs, 0.0);

    _multiply(norbs, norbs, nao, transposed.data(), dense_w.data(), square.data());

    auto closed = CPackedMatrix(norbs, norbs, mat_t::symmetric);

    for (size_t i = 0; i < norbs; i++)
    {
        for (size_t j = 0; j <= i; j++)
        {
            closed.data()[closed.index(i, j)] = square[i * norbs + j];
        }
    }

    return closed;
}

auto
CSimdRIJKGradientDriver::_check_metric(const CPackedMatrix &metric, const CMolecularBasis &aux_basis) const -> void
{
    // NOTE: the gradient needs the inverse of the metric applied twice, once as
    // the factor and once as its transpose. **Either inversion gives that**: M
    // transposed times M is the inverse of the metric for the inverted Cholesky
    // factor and for the inverted square root alike, which is why the energy takes
    // either. This driver once refused the root, and said it was because the two do
    // not close the same sum; they do. What the refusal was really protecting was a
    // multiply which skipped everything below the diagonal, and that now reads the
    // type and takes the whole matrix where it has to.
    //
    // Which of them it was handed is read from the type -- a factor is lower
    // triangular and a root is symmetric -- rather than from a flag beside it, the
    // way the Fock driver reads it.
    //
    // A range separated gradient has no choice in the matter: the attenuated metric
    // is numerically singular and its factorization fails on some fitting sets and
    // succeeds on others, so both kinds arrive.

    errors::assertMsgCritical(
        metric.number_of_rows() > 0,
        std::string("SimdRIJKGradientDriver: The metric is empty. The direct mode forms no inverted "
                    "factor, and this driver needs the one the mode which holds the B vectors forms"));

    errors::assertMsgCritical(
        metric.number_of_rows() == aux_basis.dimensions_of_basis(),
        std::string("SimdRIJKGradientDriver: The metric is not of the auxiliary basis"));
}

auto
CSimdRIJKGradientDriver::_gram(const std::vector<CPackedMatrix> &orbital_densities,
                               const size_t                      naux,
                               const size_t                      norbs) const -> std::vector<double>
{
    auto gram = std::vector<double>();

        // the elements of one packed matrix of the orbitals, which is what the
        // Gram product sums over

        const auto nelements = orbital_densities.front().number_of_elements();

        gram.assign(naux * naux, 0.0);

        auto weight = std::vector<double>(nelements, 0.0);

        for (size_t i = 0; i < norbs; i++)
        {
            for (size_t j = 0; j <= i; j++)
            {
                weight[i * (i + 1) / 2 + j] = (i == j) ? 1.0 : std::sqrt(2.0);
            }
        }

        // NOTE: the elements in panels, as the factor above is, so the weighted
        // copy is bounded by the budget. The result it adds into is the square of
        // the auxiliary basis and cannot be divided, so it is taken off first.

        const auto fixed = naux * naux * sizeof(double);

        const auto spare = (_budget > fixed) ? _budget - fixed : size_t{0};

        const auto per_element = naux * sizeof(double);

        const auto by_memory = spare / std::max(per_element, size_t{1});

        const auto npanel = std::min(nelements, std::max(_min_batch, by_memory));

        auto scaled = std::vector<double>(naux * npanel, 0.0);

        const auto nrange = static_cast<int>(naux);

        for (size_t first = 0; first < nelements; first += npanel)
        {
            const auto count = std::min(npanel, nelements - first);

#pragma omp parallel for schedule(static)
            for (int at = 0; at < nrange; at++)
            {
                const auto q = static_cast<size_t>(at);

                const auto *values = orbital_densities[q].data() + first;

                for (size_t c = 0; c < count; c++)
                {
                    scaled[q * count + c] = weight[first + c] * values[c];
                }
            }

            _add_multiply_by_transpose(naux, naux, count, scaled.data(), scaled.data(), gram.data());
        }

    return gram;
}

auto
CSimdRIJKGradientDriver::_orbital_densities(const CSparseTensor   &bq_vectors,
                                            const CMolecularBasis &basis,
                                            const CMolecularBasis &aux_basis,
                                            const CPackedMatrix   &metric,
                                            const CPackedMatrix   &coefficients) const
    -> std::vector<CPackedMatrix>
{
    const auto naux = aux_basis.dimensions_of_basis();

    const auto norbs = coefficients.number_of_columns();

    const auto nao = coefficients.number_of_rows();

    // NOTE: a spin which occupies nothing has no fitted densities and no exchange.
    // The hydrogen atom is one: one alpha orbital and no beta. Left to run, the
    // transformation asks the library for a product of no columns and it refuses.
    if (norbs == 0) return {};

    auto orbital_densities = std::vector<CPackedMatrix>();

    orbital_densities.reserve(naux);

    for (size_t q = 0; q < naux; q++)
    {
        orbital_densities.push_back(CPackedMatrix(norbs, norbs, mat_t::symmetric));
    }

    // NOTE: the half transformed W matrices are the basis functions times the
    // orbitals for every auxiliary function, which is the largest thing this
    // phase could hold and the one thing it must not hold whole: at two thousand
    // functions, two hundred orbitals and eight thousand auxiliary functions it
    // is twenty six gigabytes, where the fitted densities it closes them into are
    // one and a half. They are formed for a batch of the auxiliary basis, closed,
    // and the storage reused for the next batch.

    // NOTE: the transpose of the coefficients, formed once for the whole phase.
    // It is the basis functions times the orbitals, which is nothing beside the
    // matrices it multiplies, and it saves a transposed product per auxiliary
    // function and the reading of a transposed array in the inner loop.

    auto dense_c = std::vector<double>(nao * norbs, 0.0);

    coefficients.to_dense(dense_c.data());

    auto transposed = std::vector<double>(norbs * nao, 0.0);

    for (size_t mu = 0; mu < nao; mu++)
    {
        for (size_t i = 0; i < norbs; i++)
        {
            transposed[i * nao + mu] = dense_c[mu * norbs + i];
        }
    }

    const auto per_function = nao * norbs * sizeof(double);

    const auto by_memory = std::max(size_t{1}, _budget / std::max(per_function, size_t{1}));

    const auto nbatch = std::min(naux, std::max(_min_batch, by_memory));

    auto half = std::vector<CPackedMatrix>();

    for (size_t at = 0; at < nbatch; at++)
    {
        half.push_back(CPackedMatrix(nao, norbs, mat_t::general));
    }

    for (size_t first = 0; first < naux; first += nbatch)
    {
        const auto last = std::min(first + nbatch, naux);

        const auto count = last - first;

        auto functions = std::vector<size_t>(count);

        std::iota(functions.begin(), functions.end(), first);

        // NOTE: the last batch is shorter than the others and the transformation
        // fills what it is handed, so it is handed the front of the storage.

        if (count == nbatch)
        {
            _drv.compute_w_vectors(bq_vectors, basis, aux_basis, coefficients, functions, half);
        }
        else
        {
            auto tail = std::vector<CPackedMatrix>(half.begin(), half.begin() + static_cast<long>(count));

            _drv.compute_w_vectors(bq_vectors, basis, aux_basis, coefficients, functions, tail);

            std::copy(tail.begin(), tail.end(), half.begin());
        }

        const auto nrange = static_cast<int>(count);

#pragma omp parallel for schedule(static) if (nrange > 1)
        for (int at = 0; at < nrange; at++)
        {
            const auto q = static_cast<size_t>(at);

            orbital_densities[first + q] = _close_orbitals(transposed, nao, norbs, half[q]);
        }
    }

    _apply_transposed_factor(metric, orbital_densities);

    return orbital_densities;
}

// NOTE: the distributed phases. What divides is the auxiliary index: a rank owns
// the functions of the atoms whose B vectors it holds, forms their fitted densities
// from its own integrals, and never sees another rank's. The one place the index
// does not divide is the transposed factor of the metric, which reaches across the
// whole of it, and that is the one place the ranks talk. See the notes on
// TDistributedFit and on mpi_panel_partial.

auto
CSimdRIJKGradientDriver::_local_densities_for(const CSparseTensor        &bq_vectors,
                                              const CMolecularBasis      &basis,
                                              const CMolecularBasis      &aux_basis,
                                              const CPackedMatrix        &coefficients,
                                              const std::vector<size_t>  &functions,
                                              const size_t                naux,
                                              const size_t                budget,
                                              std::vector<CPackedMatrix> &target) const -> void
{
    const auto norbs = coefficients.number_of_columns();

    const auto nao = coefficients.number_of_rows();

    target.clear();

    // NOTE: a spin which occupies nothing has no fitted densities and no exchange.

    if ((norbs == 0) || functions.empty()) return;

    target.resize(naux);

    auto transposed = std::vector<double>(norbs * nao, 0.0);

    {
        auto dense_c = std::vector<double>(nao * norbs, 0.0);

        coefficients.to_dense(dense_c.data());

        for (size_t mu = 0; mu < nao; mu++)
        {
            for (size_t ii = 0; ii < norbs; ii++) transposed[ii * nao + mu] = dense_c[mu * norbs + ii];
        }
    }

    const auto per_function = nao * norbs * sizeof(double);

    const auto batch_by_memory = std::max(size_t{1}, budget / std::max(per_function, size_t{1}));

    const auto nbatch = std::min(functions.size(), std::max(_min_batch, batch_by_memory));

    auto half = std::vector<CPackedMatrix>();

    for (size_t at = 0; at < nbatch; at++)
    {
        half.push_back(CPackedMatrix(nao, norbs, mat_t::general));
    }

    for (size_t first = 0; first < functions.size(); first += nbatch)
    {
        const auto last = std::min(first + nbatch, functions.size());

        const auto count = last - first;

        auto wanted = std::vector<size_t>(functions.begin() + static_cast<long>(first),
                                          functions.begin() + static_cast<long>(last));

        if (count == nbatch)
        {
            _drv.compute_w_vectors(bq_vectors, basis, aux_basis, coefficients, wanted, half);
        }
        else
        {
            auto tail = std::vector<CPackedMatrix>(half.begin(), half.begin() + static_cast<long>(count));

            _drv.compute_w_vectors(bq_vectors, basis, aux_basis, coefficients, wanted, tail);

            std::copy(tail.begin(), tail.end(), half.begin());
        }

        const auto nrange = static_cast<int>(count);

#pragma omp parallel for schedule(static) if (nrange > 1)
        for (int at = 0; at < nrange; at++)
        {
            const auto local = static_cast<size_t>(at);

            target[functions[first + local]] = _close_orbitals(transposed, nao, norbs, half[local]);
        }
    }
}

auto
CSimdRIJKGradientDriver::mpi_local_densities(const CMolecule        &molecule,
                                             const CMolecularBasis  &basis,
                                             const CMolecularBasis  &aux_basis,
                                             const CSparseTensor    &bq_vectors,
                                             const CPackedMatrix    &density,
                                             const CPackedMatrix    &coefficients,
                                             const CPackedMatrix    &coefficients_beta,
                                             const std::vector<int> &aux_atoms,
                                             const size_t            budget) const -> TDistributedFit
{
    auto fit = TDistributedFit();

    fit.naux = aux_basis.dimensions_of_basis();

    const auto norbs = coefficients.number_of_columns();

    // NOTE: **empty means empty here.** aux_functions_of reads a list of no atoms
    // as every atom, which is right for a caller asking about a whole molecule and
    // wrong for every caller of this routine: these are the atoms a rank was dealt,
    // and a rank dealt none owns no functions. Expanded instead, such a rank asked
    // the transformation for every auxiliary function in the molecule while holding
    // the B vectors of none of them.
    //
    // NOTE: the third place this same convention has bitten -- the Coulomb only
    // gradient driver, the fitted Fock driver, and here.

    fit.functions = aux_atoms.empty() ? std::vector<size_t>()
                                      : simdri::aux_functions_of(aux_basis, aux_atoms);

    // the right hand side of the fitting, summed over this rank's B vectors alone.
    // The caller adds the ranks' and hands the sum back to mpi_set_fitting.

    fit.fitting = _drv.compute_y_vector(bq_vectors, basis, aux_basis, density);

    if (norbs == 0) return fit;

    fit.nelements = norbs * (norbs + 1) / 2;

    // NOTE: the panel is what the ranks exchange and what a rank holds of the whole
    // auxiliary basis at once, so it is sized from the budget this phase was given
    // rather than from the driver's. Under a communicator the ranks share a node and
    // the budget of one of them is not the budget of the machine.

    const auto per_element = 2 * fit.naux * sizeof(double);

    const auto by_memory = budget / std::max(per_element, size_t{1});

    fit.npanel = std::min(fit.nelements, std::max(_min_batch, by_memory));

    // NOTE: **both spins' panel geometry is settled before anything returns.** The
    // second spin's used to be worked out further down, after the early return
    // below, so a rank which owned no auxiliary functions came back reporting no
    // beta panels while every other rank reported some. The caller loops over the
    // panels and reduces each one, so those ranks ran different numbers of
    // collectives and the job sat in Allreduce until it was killed -- which reads
    // as water taking minutes, not as the deadlock it is.

    if (coefficients_beta.number_of_columns() > 0)
    {
        const auto norbs_beta = coefficients_beta.number_of_columns();

        fit.nelements_beta = norbs_beta * (norbs_beta + 1) / 2;

        fit.npanel_beta = std::min(fit.nelements_beta, std::max(_min_batch, by_memory));
    }

    fit.densities.resize(fit.naux);

    fit.gram_rows.assign(fit.functions.size() * fit.naux, 0.0);

    if (fit.functions.empty()) return fit;

    // NOTE: the same B vectors serve both spins. What differs between them is the
    // orbitals the half transformed vectors are closed into, so a second spin is a
    // second closing and not a second pass over the integrals.

    _local_densities_for(bq_vectors, basis, aux_basis, coefficients, fit.functions, fit.naux, budget,
                         fit.densities);

    if (coefficients_beta.number_of_columns() > 0)
    {
        errors::assertMsgCritical(coefficients_beta.number_of_elements() == 0 ||
                                      (coefficients_beta.number_of_rows() == coefficients.number_of_rows()),
                                  std::string("SimdRIJKGradientDriver: The two spins' orbitals are not of "
                                              "the same basis"));

        _local_densities_for(bq_vectors, basis, aux_basis, coefficients_beta, fit.functions, fit.naux, budget,
                             fit.densities_beta);
    }

    return fit;
}

auto
CSimdRIJKGradientDriver::mpi_set_fitting(TDistributedFit           &fit,
                                         const CPackedMatrix       &metric,
                                         const std::vector<double> &total) const -> void
{
    errors::assertMsgCritical(total.size() == fit.naux,
                              std::string("SimdRIJKGradientDriver: The summed fitting is not of the "
                                          "auxiliary basis"));

    fit.fitting = total;

    _apply_transposed_factor(metric, fit.fitting.data(), 1);
}

auto
CSimdRIJKGradientDriver::mpi_clear_fitting(TDistributedFit &fit) const -> void
{
    fit.fitting.clear();
}

auto
CSimdRIJKGradientDriver::_transposed_columns(const CPackedMatrix       &metric,
                                             const std::vector<size_t> &columns,
                                             const size_t               naux) const -> std::vector<double>
{
    // NOTE: the columns of the transposed factor this rank owns, and not the whole
    // of it. The serial phase writes out the square of the auxiliary basis; here
    // only the columns of the owned functions are ever multiplied, and on eight
    // ranks that is an eighth of the array.

    const auto triangular = (metric.get_type() != mat_t::symmetric);

    auto factor = std::vector<double>(naux * naux, 0.0);

    metric.to_dense(factor.data());

    auto columns_of = std::vector<double>(naux * columns.size(), 0.0);

    const auto nrange = static_cast<int>(naux);

#pragma omp parallel for schedule(static)
    for (int at = 0; at < nrange; at++)
    {
        const auto p = static_cast<size_t>(at);

        for (size_t k = 0; k < columns.size(); k++)
        {
            const auto q = columns[k];

            // the transpose: the element at (p, q) of it is the element at (q, p) of
            // the factor, and a triangular factor has nothing below the diagonal.

            if (triangular && (q < p)) continue;

            columns_of[p * columns.size() + k] = factor[q * naux + p];
        }
    }

    return columns_of;
}

auto
CSimdRIJKGradientDriver::mpi_panel_partial(const TDistributedFit &fit,
                                           const CPackedMatrix   &metric,
                                           const size_t           ipanel,
                                           const bool             beta) const -> std::vector<double>
{
    const auto range = fit.panel_range(ipanel, beta);

    const auto first = range.first;

    const auto count = range.second;

    if (count == 0) return {};

    auto partial = std::vector<double>(fit.naux * count, 0.0);

    const auto &densities = beta ? fit.densities_beta : fit.densities;

    if (fit.functions.empty() || densities.empty()) return partial;

    // this rank's rows of the panel, gathered so that the owned auxiliary index is
    // the one a single product contracts

    auto gathered = std::vector<double>(fit.functions.size() * count, 0.0);

    const auto nowned = static_cast<int>(fit.functions.size());

#pragma omp parallel for schedule(static)
    for (int at = 0; at < nowned; at++)
    {
        const auto k = static_cast<size_t>(at);

        const auto &matrix = densities[fit.functions[k]];

        if (matrix.number_of_elements() == 0) continue;

        std::copy_n(matrix.data() + first, count, gathered.data() + k * count);
    }

    const auto columns = _transposed_columns(metric, fit.functions, fit.naux);

    _multiply(fit.naux, count, fit.functions.size(), columns.data(), gathered.data(), partial.data());

    return partial;
}

auto
CSimdRIJKGradientDriver::mpi_panel_absorb(TDistributedFit &fit,
                                          const size_t     ipanel,
                                          const double    *reduced,
                                          const size_t     size,
                                          const bool       beta) const -> void
{
    const auto range = fit.panel_range(ipanel, beta);

    const auto first = range.first;

    const auto count = range.second;

    if (count == 0) return;

    errors::assertMsgCritical(size == fit.naux * count,
                              std::string("SimdRIJKGradientDriver: The reduced panel is not the auxiliary "
                                          "basis by the elements of the panel"));

    auto &densities = beta ? fit.densities_beta : fit.densities;

    // the rows this rank owns, written back over the ones it formed. What was there
    // was the fitted density before the metric; what goes in is after it.

    const auto nowned = static_cast<int>(fit.functions.size());

#pragma omp parallel for schedule(static)
    for (int at = 0; at < nowned; at++)
    {
        const auto k = static_cast<size_t>(at);

        auto &matrix = densities[fit.functions[k]];

        if (matrix.number_of_elements() == 0) continue;

        std::copy_n(reduced + fit.functions[k] * count, count, matrix.data() + first);
    }

    // and this rank's rows of the Gram, from the whole of the panel while it is
    // here. The weights are the ones the serial phase uses: an element off the
    // diagonal of a packed matrix stands for two of the sum over the pairs of
    // orbitals, and the root of that multiplies once from each side.

    // NOTE: the two spins add into one Gram, as they do in the serial phase, so the
    // second is accumulated into the same rows rather than into rows of its own.

    // NOTE: the weights are written out over the whole of the elements and the
    // panel's slice taken, which is the same double loop the serial Gram uses. It
    // is the orbitals squared over two doubles and is nothing beside the products
    // below; working the row of an element out from the element instead would be an
    // inverse triangular number, and getting it wrong would weight the Gram wrongly
    // in a way no dimension check would catch.

    size_t norbs = 0;

    for (const auto q : fit.functions)
    {
        if (densities[q].number_of_elements() > 0)
        {
            norbs = densities[q].number_of_rows();

            break;
        }
    }

    if (norbs == 0) return;

    auto all_weights = std::vector<double>(beta ? fit.nelements_beta : fit.nelements, 0.0);

    for (size_t i = 0; i < norbs; i++)
    {
        for (size_t j = 0; j <= i; j++)
        {
            all_weights[i * (i + 1) / 2 + j] = (i == j) ? 1.0 : std::sqrt(2.0);
        }
    }

    auto weight = std::vector<double>(all_weights.begin() + static_cast<long>(first),
                                      all_weights.begin() + static_cast<long>(first + count));

    auto scaled = std::vector<double>(fit.naux * count, 0.0);

    const auto nrange = static_cast<int>(fit.naux);

#pragma omp parallel for schedule(static)
    for (int at = 0; at < nrange; at++)
    {
        const auto q = static_cast<size_t>(at);

        for (size_t c = 0; c < count; c++)
        {
            scaled[q * count + c] = weight[c] * reduced[q * count + c];
        }
    }

    // rows of this rank against the whole: gram_rows(|owned| x naux) += owned(|owned| x count) * scaled^T

    auto owned = std::vector<double>(fit.functions.size() * count, 0.0);

#pragma omp parallel for schedule(static)
    for (int at = 0; at < nowned; at++)
    {
        const auto k = static_cast<size_t>(at);

        std::copy_n(scaled.data() + fit.functions[k] * count, count, owned.data() + k * count);
    }

    _add_multiply_by_transpose(fit.functions.size(), fit.naux, count, owned.data(), scaled.data(),
                               fit.gram_rows.data());
}

auto
CSimdRIJKGradientDriver::mpi_omega(const TDistributedFit &fit,
                                   const double          *gram,
                                   const size_t           size,
                                   const double           exchange_scaling_factor,
                                   const bool             open_shell) const -> CPackedMatrix
{
    const auto naux = fit.naux;

    auto omega = CPackedMatrix(naux, naux, mat_t::symmetric);

    omega.zero();

    const auto has_gram = (gram != nullptr) && (size > 0) && (exchange_scaling_factor != 0.0);

    errors::assertMsgCritical((!has_gram) || (size == naux * naux),
                              std::string("SimdRIJKGradientDriver: The gathered Gram is not the square of "
                                          "the auxiliary basis"));

    // NOTE: the same factors the serial phases carry, and for the same reasons. A
    // closed shell is handed one spin's density and carries two on the Coulomb and
    // one on a single spin's Gram; an open shell is handed the total and carries a
    // half on each, the two spins having been added into one Gram already.

    const auto coulomb_factor = open_shell ? 0.5 : 2.0;

    const auto exchange_factor = open_shell ? 0.5 * exchange_scaling_factor : exchange_scaling_factor;

    // NOTE: an empty fitting is how the attenuated operator says it has no Coulomb
    // term, which is the same signal the serial phase uses: it returns the fitting
    // empty rather than as zeros, so that a caller which reaches for it gets an
    // error and not a number which looks like an answer.

    const auto has_coulomb = (!fit.fitting.empty());

    for (size_t p = 0; p < naux; p++)
    {
        for (size_t q = 0; q <= p; q++)
        {
            auto value = has_coulomb ? coulomb_factor * fit.fitting[p] * fit.fitting[q] : 0.0;

            if (has_gram) value -= exchange_factor * gram[p * naux + q];

            omega.data()[omega.index(p, q)] = value;
        }
    }

    return omega;
}

auto
CSimdRIJKGradientDriver::mpi_compute_share(const CMolecule        &molecule,
                                           const CMolecularBasis  &basis,
                                           const CMolecularBasis  &aux_basis,
                                           const TDistributedFit  &fit,
                                           const CPackedMatrix    &density,
                                           const CPackedMatrix    &coefficients,
                                           const CPackedMatrix    &coefficients_beta,
                                           const CPackedMatrix    &omega,
                                           const double            exchange_scaling_factor,
                                           const std::vector<int> &atoms,
                                           const std::vector<int> &aux_atoms) const -> CPackedMatrix
{
    const auto natoms = molecule.number_of_atoms();

    auto gradient = CPackedMatrix(natoms, 3, mat_t::general);

    gradient.zero();

    auto wanted = std::vector<bool>(natoms, false);

    for (const auto iatom : atoms)
    {
        errors::assertMsgCritical((iatom >= 0) && (static_cast<size_t>(iatom) < natoms),
                                  std::string("SimdRIJKGradientDriver: An atom outside the molecule was "
                                              "asked for"));

        wanted[static_cast<size_t>(iatom)] = true;
    }

    // NOTE: the fitted densities are the length of the whole auxiliary basis with
    // only this rank's entries filled, and the contraction below reads them only for
    // the auxiliary atoms it is given, which are this rank's. An entry of another
    // rank's is never touched, which is what lets the serial contraction stand
    // unchanged under a communicator.

    const auto open_shell = (!fit.densities_beta.empty());

    // NOTE: each spin with its own orbitals. Handing the first spin's for both is
    // the kind of mistake which leaves a gradient that looks like a gradient.

    const auto spins = open_shell
                           ? std::vector<TExchangeSpin>{{&coefficients, &fit.densities},
                                                        {&coefficients_beta, &fit.densities_beta}}
                           : std::vector<TExchangeSpin>{{&coefficients, &fit.densities}};

    const auto coulomb_factor = open_shell ? 1.0 : 4.0;

    const auto exchange_factor = open_shell ? -exchange_scaling_factor : -2.0 * exchange_scaling_factor;

    // NOTE: a rank dealt no auxiliary atoms adds no three-center term. It is not
    // passed on to the contraction, which reads an empty list as every atom.

    if (!aux_atoms.empty())
    {
        _compute_three_center(gradient, molecule, basis, aux_basis, fit.fitting, density, spins,
                              coulomb_factor, exchange_factor, wanted, aux_atoms);
    }

    // NOTE: the two-center term is asked of one rank alone. Omega is of the whole
    // auxiliary basis on every rank, so every rank could form it and the ranks would
    // then add it as many times as there are of them. An empty matrix is how a rank
    // says it is not the one.

    if (omega.number_of_elements() > 0)
    {
        const auto two_center = CSimdTwoCenterElectronRepulsionGradientDriver(_block_size);

        const auto metric_part = two_center.compute(molecule, aux_basis, omega, atoms);

        for (size_t iatom = 0; iatom < natoms; iatom++)
        {
            for (size_t c = 0; c < 3; c++)
            {
                gradient.data()[gradient.index(iatom, c)] -= metric_part.at(iatom, c);
            }
        }
    }

    return gradient;
}

auto
CSimdRIJKGradientDriver::mpi_compute_share_rs(const CMolecule        &molecule,
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
                                              const std::vector<int> &aux_atoms) const -> CPackedMatrix
{
    errors::assertMsgCritical(omega > 0.0,
                              std::string("SimdRIJKGradientDriver: A range separated share was asked for "
                                          "without a range separation parameter"));

    const auto natoms = molecule.number_of_atoms();

    auto gradient = CPackedMatrix(natoms, 3, mat_t::general);

    gradient.zero();

    auto wanted = std::vector<bool>(natoms, false);

    for (const auto iatom : atoms)
    {
        errors::assertMsgCritical((iatom >= 0) && (static_cast<size_t>(iatom) < natoms),
                                  std::string("SimdRIJKGradientDriver: An atom outside the molecule was "
                                              "asked for"));

        wanted[static_cast<size_t>(iatom)] = true;
    }

    const auto open_shell = (!fit.densities_beta.empty());

    const auto spins = open_shell
                           ? std::vector<TExchangeSpin>{{&coefficients, &fit.densities},
                                                        {&coefficients_beta, &fit.densities_beta}}
                           : std::vector<TExchangeSpin>{{&coefficients, &fit.densities}};

    const auto spins_erf = open_shell
                               ? std::vector<TExchangeSpin>{{&coefficients, &fit_erf.densities},
                                                            {&coefficients_beta, &fit_erf.densities_beta}}
                               : std::vector<TExchangeSpin>{{&coefficients, &fit_erf.densities}};

    // NOTE: the factors of the serial entry. The two operators come out of one call
    // of the derivative driver, so the attenuated term costs a second contraction
    // and not a second pass over the integrals; and the attenuated operator has no
    // Coulomb term, which is why only one Coulomb factor appears.

    const auto coulomb_factor = open_shell ? 1.0 : 4.0;

    const auto exchange_factor = open_shell ? -exchange_scaling_factor : -2.0 * exchange_scaling_factor;

    const auto erf_factor = open_shell ? -erf_exchange_scaling_factor : -2.0 * erf_exchange_scaling_factor;

    if (!aux_atoms.empty())
    {
        _compute_three_center(gradient, molecule, basis, aux_basis, fit.fitting, density, spins,
                              coulomb_factor, exchange_factor, wanted, aux_atoms, spins_erf, erf_factor,
                              omega);
    }

    if (omega_plain.number_of_elements() > 0)
    {
        const auto two_center = CSimdTwoCenterElectronRepulsionGradientRsDriver(_block_size);

        // NOTE: the attenuated driver answers with the two terms apart, the plain
        // operator's and the attenuated one's, and both are subtracted.

        const auto parts = two_center.compute(molecule, aux_basis, omega_plain, omega_erf, omega, atoms);

        for (size_t iatom = 0; iatom < natoms; iatom++)
        {
            for (size_t c = 0; c < 3; c++)
            {
                gradient.data()[gradient.index(iatom, c)] -= parts.first.at(iatom, c);

                gradient.data()[gradient.index(iatom, c)] -= parts.second.at(iatom, c);
            }
        }
    }

    return gradient;
}

auto
CSimdRIJKGradientDriver::fitted_densities(const CSparseTensor   &bq_vectors,
                                          const CMolecularBasis &basis,
                                          const CMolecularBasis &aux_basis,
                                          const CPackedMatrix   &metric,
                                          const CPackedMatrix   &density,
                                          const CPackedMatrix   &coefficients,
                                          const double           exchange_scaling_factor) const
    -> TFittedDensities
{
    _check_metric(metric, aux_basis);

    const auto naux = aux_basis.dimensions_of_basis();

    const auto norbs = coefficients.number_of_columns();

    // the fitting coefficients. compute_y_vector closes the sum over the basis
    // functions of the B vectors against the density, which is the factor already
    // applied once; the transpose of it applied to that is the inverse.

    auto fitting = _drv.compute_y_vector(bq_vectors, basis, aux_basis, density);

    _apply_transposed_factor(metric, fitting.data(), 1);

    // the fitted densities of the occupied orbitals. The first index is
    // transformed by the driver which holds the B vectors and the second here,
    // which leaves a symmetric matrix of the orbitals for each auxiliary
    // function; the metric is applied to those and never in the basis of the
    // atomic orbitals.

    auto orbital_densities = _orbital_densities(bq_vectors, basis, aux_basis, metric, coefficients);

    // and the two index fitted density, which the derivative of the metric is
    // contracted against. The exchange part is the Gram product of the fitted
    // densities of the orbitals and the Coulomb part is the outer product of the
    // fitting coefficients.

    auto omega = CPackedMatrix(naux, naux, mat_t::symmetric);

    omega.zero();

    // NOTE: the exchange part is the Gram product of the fitted densities over
    // every pair of auxiliary functions, which as a sum written out is the square
    // of the auxiliary basis times the square of the orbitals and was the whole
    // of this phase on a molecule of any size, on one thread. It is a product of
    // one matrix by its own transpose, so it is taken as one.
    //
    // The packed matrices hold the lower triangle, so an element off the diagonal
    // stands for two of the sum over the pairs of orbitals. Scaling each element
    // by the root of what it counts for puts that weight into the product: the
    // root multiplies twice, once from each side, and returns the two.

    auto gram = std::vector<double>();

    if (exchange_scaling_factor != 0.0)
    {
        gram = _gram(orbital_densities, naux, norbs);
    }

    for (size_t p = 0; p < naux; p++)
    {
        for (size_t q = 0; q <= p; q++)
        {
            auto value = 2.0 * fitting[p] * fitting[q];

            if (!gram.empty())
            {
                value -= exchange_scaling_factor * gram[p * naux + q];
            }

            omega.data()[omega.index(p, q)] = value;
        }
    }

    return {std::move(fitting), std::move(orbital_densities), {}, std::move(omega)};
}

auto
CSimdRIJKGradientDriver::compute_rs(const CMolecule        &molecule,
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
                                    const std::vector<int> &aux_atoms) const -> CPackedMatrix
{
    const auto natoms = molecule.number_of_atoms();

    for (const auto iatom : atoms)
    {
        errors::assertMsgCritical(
            (iatom >= 0) && (static_cast<size_t>(iatom) < natoms),
            std::string("SimdRIJKGradientDriver.compute_rs: Atom is not an atom of the molecule"));
    }

    errors::assertMsgCritical(
        std::set<int>(atoms.begin(), atoms.end()).size() == atoms.size(),
        std::string("SimdRIJKGradientDriver.compute_rs: An atom is named more than once"));

    errors::assertMsgCritical(
        density.get_type() == mat_t::symmetric,
        std::string("SimdRIJKGradientDriver.compute_rs: The density matrix is expected to be symmetric"));

    errors::assertMsgCritical(
        omega > 0.0,
        std::string("SimdRIJKGradientDriver.compute_rs: The range separation parameter must be positive. A "
                    "functional which is not range separated is served by compute"));

    auto gradient = CPackedMatrix(natoms, 3, mat_t::general);

    gradient.zero();

    if (atoms.empty()) return gradient;

    auto wanted = std::vector<bool>(natoms, false);

    for (const auto iatom : atoms) wanted[static_cast<size_t>(iatom)] = true;

    // NOTE: two sets of fitted densities. The plain one carries the Coulomb term
    // and the plain exchange; the attenuated one carries its exchange alone, is
    // fitted in its own metric, and comes back with an empty fitting.

    const auto fitted = fitted_densities(bq_vectors, basis, aux_basis, metric, density, coefficients,
                                         exchange_scaling_factor);

    const auto fitted_erf = fitted_densities_rs(bq_vectors_erf, basis, aux_basis, metric_erf, coefficients,
                                                erf_exchange_scaling_factor);

    // NOTE: the factors of a closed shell, as in compute above: the density handed
    // in is one spin's, so the Coulomb carries four and each exchange twice its
    // coefficient.

    const auto spins = std::vector<TExchangeSpin>{{&coefficients, &fitted.orbital_densities}};

    const auto spins_erf = std::vector<TExchangeSpin>{{&coefficients, &fitted_erf.orbital_densities}};

    _compute_three_center(gradient, molecule, basis, aux_basis, fitted.coefficients, density, spins,
                          4.0, -2.0 * exchange_scaling_factor, wanted, aux_atoms, spins_erf,
                          -2.0 * erf_exchange_scaling_factor, omega);

    // NOTE: both metric terms, from one call of the two-center range separated
    // driver. Each operator's derivative is contracted with its own weighted
    // density; crossing them would be a gradient of nothing.

    const auto two_center = CSimdTwoCenterElectronRepulsionGradientRsDriver(_block_size);

    const auto [metric_part, metric_part_erf] =
        two_center.compute(molecule, aux_basis, fitted.omega, fitted_erf.omega, omega, atoms);

    for (size_t iatom = 0; iatom < natoms; iatom++)
    {
        for (size_t c = 0; c < 3; c++)
        {
            gradient.data()[gradient.index(iatom, c)] -= metric_part.at(iatom, c);

            gradient.data()[gradient.index(iatom, c)] -= metric_part_erf.at(iatom, c);
        }
    }

    return gradient;
}

auto
CSimdRIJKGradientDriver::compute_rs(const CMolecule       &molecule,
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
                                    const double           omega) const -> CPackedMatrix
{
    auto atoms = std::vector<int>(molecule.number_of_atoms());

    std::iota(atoms.begin(), atoms.end(), 0);

    return compute_rs(molecule, basis, aux_basis, bq_vectors, bq_vectors_erf, metric, metric_erf, density,
                      coefficients, exchange_scaling_factor, erf_exchange_scaling_factor, omega, atoms, {});
}

auto
CSimdRIJKGradientDriver::compute_open_shell_rs(const CMolecule        &molecule,
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
                                               const std::vector<int> &aux_atoms) const -> CPackedMatrix
{
    const auto natoms = static_cast<size_t>(molecule.number_of_atoms());

    for (const auto iatom : atoms)
    {
        errors::assertMsgCritical(
            (iatom >= 0) && (static_cast<size_t>(iatom) < natoms),
            std::string("SimdRIJKGradientDriver.compute_open_shell_rs: Atom is not an atom of the molecule"));
    }

    errors::assertMsgCritical(
        std::set<int>(atoms.begin(), atoms.end()).size() == atoms.size(),
        std::string("SimdRIJKGradientDriver.compute_open_shell_rs: An atom is named more than once"));

    errors::assertMsgCritical(
        density.get_type() == mat_t::symmetric,
        std::string("SimdRIJKGradientDriver.compute_open_shell_rs: The density matrix is expected to be symmetric"));

    errors::assertMsgCritical(
        omega > 0.0,
        std::string("SimdRIJKGradientDriver.compute_open_shell_rs: The range separation parameter must be "
                    "positive. A functional which is not range separated is served by compute_open_shell"));

    auto gradient = CPackedMatrix(natoms, 3, mat_t::general);

    gradient.zero();

    if (atoms.empty()) return gradient;

    auto wanted = std::vector<bool>(natoms, false);

    for (const auto iatom : atoms) wanted[static_cast<size_t>(iatom)] = true;

    const auto fitted = fitted_densities_open_shell(bq_vectors, basis, aux_basis, metric, density,
                                                    coefficients_alpha, coefficients_beta,
                                                    exchange_scaling_factor);

    const auto fitted_erf = fitted_densities_open_shell_rs(bq_vectors_erf, basis, aux_basis, metric_erf,
                                                           coefficients_alpha, coefficients_beta,
                                                           erf_exchange_scaling_factor);

    // NOTE: the factors of an open shell, as in compute_open_shell above: the
    // density is the total one, so the Coulomb carries one and each spin's
    // exchange carries its coefficient once.

    const auto spins = std::vector<TExchangeSpin>{{&coefficients_alpha, &fitted.orbital_densities},
                                                  {&coefficients_beta, &fitted.orbital_densities_beta}};

    const auto spins_erf = std::vector<TExchangeSpin>{{&coefficients_alpha, &fitted_erf.orbital_densities},
                                                      {&coefficients_beta, &fitted_erf.orbital_densities_beta}};

    _compute_three_center(gradient, molecule, basis, aux_basis, fitted.coefficients, density, spins,
                          1.0, -exchange_scaling_factor, wanted, aux_atoms, spins_erf,
                          -erf_exchange_scaling_factor, omega);

    const auto two_center = CSimdTwoCenterElectronRepulsionGradientRsDriver(_block_size);

    const auto [metric_part, metric_part_erf] =
        two_center.compute(molecule, aux_basis, fitted.omega, fitted_erf.omega, omega, atoms);

    for (size_t iatom = 0; iatom < natoms; iatom++)
    {
        for (size_t c = 0; c < 3; c++)
        {
            gradient.data()[gradient.index(iatom, c)] -= metric_part.at(iatom, c);

            gradient.data()[gradient.index(iatom, c)] -= metric_part_erf.at(iatom, c);
        }
    }

    return gradient;
}

auto
CSimdRIJKGradientDriver::fitted_densities_rs(const CSparseTensor   &bq_vectors_erf,
                                             const CMolecularBasis &basis,
                                             const CMolecularBasis &aux_basis,
                                             const CPackedMatrix   &metric_erf,
                                             const CPackedMatrix   &coefficients,
                                             const double           erf_exchange_scaling_factor) const
    -> TFittedDensities
{
    _check_metric(metric_erf, aux_basis);

    const auto naux = aux_basis.dimensions_of_basis();

    const auto norbs = coefficients.number_of_columns();

    // NOTE: no fitting coefficients and no Coulomb part of Omega. The attenuated
    // operator appears in the Fock matrix only through the exchange; the Coulomb
    // term of a range separated functional is the whole of 1 / r and is fitted in
    // the plain metric by the routine above. Fitting the density in the attenuated
    // metric would be a well formed calculation of a quantity nothing wants.

    auto orbital_densities = _orbital_densities(bq_vectors_erf, basis, aux_basis, metric_erf, coefficients);

    auto omega = CPackedMatrix(naux, naux, mat_t::symmetric);

    omega.zero();

    if ((erf_exchange_scaling_factor != 0.0) && (norbs > 0))
    {
        const auto gram = _gram(orbital_densities, naux, norbs);

        for (size_t p = 0; p < naux; p++)
        {
            for (size_t q = 0; q <= p; q++)
            {
                omega.data()[omega.index(p, q)] = -erf_exchange_scaling_factor * gram[p * naux + q];
            }
        }
    }

    // NOTE: the fitting is returned empty rather than as zeros, so that a caller
    // which passes it to the three-center term with a Coulomb factor which is not
    // zero is stopped rather than quietly given nothing.

    return {std::vector<double>(), std::move(orbital_densities), {}, std::move(omega)};
}

auto
CSimdRIJKGradientDriver::fitted_densities_open_shell_rs(const CSparseTensor   &bq_vectors_erf,
                                                        const CMolecularBasis &basis,
                                                        const CMolecularBasis &aux_basis,
                                                        const CPackedMatrix   &metric_erf,
                                                        const CPackedMatrix   &coefficients_alpha,
                                                        const CPackedMatrix   &coefficients_beta,
                                                        const double           erf_exchange_scaling_factor) const
    -> TFittedDensities
{
    _check_metric(metric_erf, aux_basis);

    const auto naux = aux_basis.dimensions_of_basis();

    const auto norbs_alpha = coefficients_alpha.number_of_columns();

    const auto norbs_beta = coefficients_beta.number_of_columns();

    auto alpha = _orbital_densities(bq_vectors_erf, basis, aux_basis, metric_erf, coefficients_alpha);

    auto beta = _orbital_densities(bq_vectors_erf, basis, aux_basis, metric_erf, coefficients_beta);

    auto omega = CPackedMatrix(naux, naux, mat_t::symmetric);

    omega.zero();

    if (erf_exchange_scaling_factor != 0.0)
    {
        auto gram = std::vector<double>();

        if (norbs_alpha > 0) gram = _gram(alpha, naux, norbs_alpha);

        if (norbs_beta > 0)
        {
            auto other = _gram(beta, naux, norbs_beta);

            if (gram.empty())
            {
                gram = std::move(other);
            }
            else
            {
                for (size_t i = 0; i < gram.size(); i++) gram[i] += other[i];
            }
        }

        // NOTE: the half which the open shell carries, as in the routine above:
        // writing the closed shell expression in terms of a sum over the spins puts
        // a half on the exchange. Setting the two spins equal returns the closed
        // shell Omega of the attenuated operator exactly, which is the first check
        // to make on this.

        if (!gram.empty())
        {
            for (size_t p = 0; p < naux; p++)
            {
                for (size_t q = 0; q <= p; q++)
                {
                    omega.data()[omega.index(p, q)] =
                        -0.5 * erf_exchange_scaling_factor * gram[p * naux + q];
                }
            }
        }
    }

    return {std::vector<double>(), std::move(alpha), std::move(beta), std::move(omega)};
}

auto
CSimdRIJKGradientDriver::fitted_densities_open_shell(const CSparseTensor   &bq_vectors,
                                                     const CMolecularBasis &basis,
                                                     const CMolecularBasis &aux_basis,
                                                     const CPackedMatrix   &metric,
                                                     const CPackedMatrix   &density,
                                                     const CPackedMatrix   &coefficients_alpha,
                                                     const CPackedMatrix   &coefficients_beta,
                                                     const double           exchange_scaling_factor) const
    -> TFittedDensities
{
    _check_metric(metric, aux_basis);

    const auto naux = aux_basis.dimensions_of_basis();

    // NOTE: the density handed in is the total one, of both spins, so the fitting
    // coefficients are of the total density. The closed shell routine is handed
    // one spin's and carries the factors of two which follow from that; here they
    // are absent and the factors below are the plain ones.

    auto fitting = _drv.compute_y_vector(bq_vectors, basis, aux_basis, density);

    _apply_transposed_factor(metric, fitting.data(), 1);

    // the fitted densities of each spin's occupied orbitals

    auto alpha = _orbital_densities(bq_vectors, basis, aux_basis, metric, coefficients_alpha);

    auto beta = _orbital_densities(bq_vectors, basis, aux_basis, metric, coefficients_beta);

    auto omega = CPackedMatrix(naux, naux, mat_t::symmetric);

    omega.zero();

    // NOTE: the Gram product of each spin, added. A spin which occupies nothing
    // contributes none, which is the hydrogen atom and is not a special case here:
    // its list of fitted densities is empty and the product is skipped.

    auto gram = std::vector<double>();

    if (exchange_scaling_factor != 0.0)
    {
        const auto norbs_alpha = coefficients_alpha.number_of_columns();

        const auto norbs_beta = coefficients_beta.number_of_columns();

        if (norbs_alpha > 0) gram = _gram(alpha, naux, norbs_alpha);

        if (norbs_beta > 0)
        {
            auto other = _gram(beta, naux, norbs_beta);

            if (gram.empty())
            {
                gram = std::move(other);
            }
            else
            {
                for (size_t i = 0; i < gram.size(); i++) gram[i] += other[i];
            }
        }
    }

    // NOTE: the two centre term of an open shell. Writing the closed shell one in
    // terms of the total density and a sum over the spins gives a half on the
    // Coulomb, where the closed shell has two against a one spin density, and a
    // half on the exchange, where it has one against a single spin's Gram. Setting
    // the two spins equal returns the closed shell expression exactly, which is
    // what the first of the checks on this exercises.

    for (size_t p = 0; p < naux; p++)
    {
        for (size_t q = 0; q <= p; q++)
        {
            auto value = 0.5 * fitting[p] * fitting[q];

            if (!gram.empty())
            {
                value -= 0.5 * exchange_scaling_factor * gram[p * naux + q];
            }

            omega.data()[omega.index(p, q)] = value;
        }
    }

    return {std::move(fitting), std::move(alpha), std::move(beta), std::move(omega)};
}

auto
CSimdRIJKGradientDriver::compute_open_shell(const CMolecule        &molecule,
                                            const CMolecularBasis  &basis,
                                            const CMolecularBasis  &aux_basis,
                                            const CSparseTensor    &bq_vectors,
                                            const CPackedMatrix    &metric,
                                            const CPackedMatrix    &density,
                                            const CPackedMatrix    &coefficients_alpha,
                                            const CPackedMatrix    &coefficients_beta,
                                            const double            exchange_scaling_factor,
                                            const std::vector<int> &atoms,
                                            const std::vector<int> &aux_atoms) const -> CPackedMatrix
{
    const auto natoms = static_cast<size_t>(molecule.number_of_atoms());

    for (const auto iatom : atoms)
    {
        errors::assertMsgCritical(
            (iatom >= 0) && (static_cast<size_t>(iatom) < natoms),
            std::string("SimdRIJKGradientDriver.compute_open_shell: Atom is not an atom of the molecule"));
    }

    errors::assertMsgCritical(
        std::set<int>(atoms.begin(), atoms.end()).size() == atoms.size(),
        std::string("SimdRIJKGradientDriver.compute_open_shell: An atom is named more than once"));

    errors::assertMsgCritical(
        density.get_type() == mat_t::symmetric,
        std::string("SimdRIJKGradientDriver.compute_open_shell: The density matrix is expected to be symmetric"));

    auto gradient = CPackedMatrix(natoms, 3, mat_t::general);

    gradient.zero();

    if (atoms.empty()) return gradient;

    auto wanted = std::vector<bool>(natoms, false);

    for (const auto iatom : atoms) wanted[static_cast<size_t>(iatom)] = true;

    const auto fitted = fitted_densities_open_shell(bq_vectors, basis, aux_basis, metric, density,
                                                    coefficients_alpha, coefficients_beta,
                                                    exchange_scaling_factor);

    // NOTE: two spins, and the factors of an open shell. The density is the total
    // one and the fitting coefficients are of it, so the Coulomb carries one where
    // the closed shell carries four, and each spin's exchange carries the fraction
    // of exact exchange once where the closed shell carries it twice for both.

    const auto spins = std::vector<TExchangeSpin>{{&coefficients_alpha, &fitted.orbital_densities},
                                                  {&coefficients_beta, &fitted.orbital_densities_beta}};

    _compute_three_center(gradient, molecule, basis, aux_basis, fitted.coefficients, density, spins,
                          1.0, -exchange_scaling_factor, wanted, aux_atoms);

    const auto two_center = CSimdTwoCenterElectronRepulsionGradientDriver(_block_size);

    const auto metric_part = two_center.compute(molecule, aux_basis, fitted.omega, atoms);

    for (size_t iatom = 0; iatom < natoms; iatom++)
    {
        for (size_t c = 0; c < 3; c++)
        {
            gradient.data()[gradient.index(iatom, c)] -= metric_part.at(iatom, c);
        }
    }

    return gradient;
}

auto
CSimdRIJKGradientDriver::compute_open_shell(const CMolecule       &molecule,
                                            const CMolecularBasis &basis,
                                            const CMolecularBasis &aux_basis,
                                            const CSparseTensor   &bq_vectors,
                                            const CPackedMatrix   &metric,
                                            const CPackedMatrix   &density,
                                            const CPackedMatrix   &coefficients_alpha,
                                            const CPackedMatrix   &coefficients_beta,
                                            const double           exchange_scaling_factor) const -> CPackedMatrix
{
    auto atoms = std::vector<int>(molecule.number_of_atoms());

    std::iota(atoms.begin(), atoms.end(), 0);

    return compute_open_shell(molecule, basis, aux_basis, bq_vectors, metric, density,
                              coefficients_alpha, coefficients_beta, exchange_scaling_factor,
                              atoms, {});
}

auto
CSimdRIJKGradientDriver::get_memory_budget() const -> size_t
{
    return _budget;
}

auto
CSimdRIJKGradientDriver::compute(const CMolecule        &molecule,
                                 const CMolecularBasis  &basis,
                                 const CMolecularBasis  &aux_basis,
                                 const CSparseTensor    &bq_vectors,
                                 const CPackedMatrix    &metric,
                                 const CPackedMatrix    &density,
                                 const CPackedMatrix    &coefficients,
                                 const double            exchange_scaling_factor,
                                 const std::vector<int> &atoms,
                                 const std::vector<int> &aux_atoms) const -> CPackedMatrix
{
    const auto natoms = molecule.number_of_atoms();

    // NOTE: an atom named twice would be accumulated twice, and an atom which is
    // not in the molecule has no row to accumulate into. Both are caught here
    // rather than by the arithmetic, which would give a gradient that is merely
    // wrong.
    for (const auto iatom : atoms)
    {
        errors::assertMsgCritical(
            (iatom >= 0) && (static_cast<size_t>(iatom) < natoms),
            std::string("SimdRIJKGradientDriver.compute: Atom is not an atom of the molecule"));
    }

    errors::assertMsgCritical(
        std::set<int>(atoms.begin(), atoms.end()).size() == atoms.size(),
        std::string("SimdRIJKGradientDriver.compute: An atom is named more than once"));

    errors::assertMsgCritical(
        density.get_type() == mat_t::symmetric,
        std::string("SimdRIJKGradientDriver.compute: The density matrix is expected to be symmetric"));

    auto gradient = CPackedMatrix(natoms, 3, mat_t::general);

    gradient.zero();

    if (atoms.empty()) return gradient;

    auto wanted = std::vector<bool>(natoms, false);

    for (const auto iatom : atoms) wanted[static_cast<size_t>(iatom)] = true;

    // the fitted densities, which read no integrals

    const auto fitted = fitted_densities(bq_vectors, basis, aux_basis, metric, density, coefficients,
                                         exchange_scaling_factor);

    // the three-center term, plus in the gradient

    // NOTE: one spin, and the factors of a closed shell: the density handed in is
    // one spin's and the fitting coefficients are of that density, so the Coulomb
    // carries four and the exchange twice the fraction of exact exchange.
    const auto spins = std::vector<TExchangeSpin>{{&coefficients, &fitted.orbital_densities}};

    _compute_three_center(gradient, molecule, basis, aux_basis, fitted.coefficients, density, spins,
                          4.0, -2.0 * exchange_scaling_factor, wanted, aux_atoms);

    // and the two-center one, which the note carries with a minus. The driver of
    // it returns the sum of Omega against the derivative with no sign applied, so
    // the sign is taken here, where it is known which term is being asked for.

    const auto two_center = CSimdTwoCenterElectronRepulsionGradientDriver(_block_size);

    const auto metric_part = two_center.compute(molecule, aux_basis, fitted.omega, atoms);

    for (size_t iatom = 0; iatom < natoms; iatom++)
    {
        for (size_t c = 0; c < 3; c++)
        {
            gradient.data()[gradient.index(iatom, c)] -= metric_part.at(iatom, c);
        }
    }

    return gradient;
}

auto
CSimdRIJKGradientDriver::compute(const CMolecule       &molecule,
                                 const CMolecularBasis &basis,
                                 const CMolecularBasis &aux_basis,
                                 const CSparseTensor   &bq_vectors,
                                 const CPackedMatrix   &metric,
                                 const CPackedMatrix   &density,
                                 const CPackedMatrix   &coefficients,
                                 const double           exchange_scaling_factor) const -> CPackedMatrix
{
    // NOTE: every atom of the molecule, which is what a gradient usually means.
    // The form which takes a list is for a caller holding a share of them.
    auto atoms = std::vector<int>(molecule.number_of_atoms());

    std::iota(atoms.begin(), atoms.end(), 0);

    return compute(molecule, basis, aux_basis, bq_vectors, metric, density, coefficients,
                   exchange_scaling_factor, atoms, {});
}

auto
CSimdRIJKGradientDriver::_compute_three_center(CPackedMatrix                    &gradient,
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
                                               const std::vector<TExchangeSpin> &spins_erf,
                                               const double                      erf_exchange_factor,
                                               const double                      omega) const -> void
{
    const auto nao = basis.dimensions_of_basis();

    const auto indices = denseidx::index_functions(basis);

    const auto aux_indices = denseidx::index_functions(aux_basis);

    const auto starts = denseidx::make_dense_starts(basis);

    const auto strides = denseidx::make_dense_strides(basis);

    const auto aux_starts = denseidx::make_dense_starts(aux_basis);

    const auto aux_strides = denseidx::make_dense_strides(aux_basis);

    const auto nmoms = strides.size();

    const auto aux_nmoms = aux_strides.size();

    // the density as a square, read once for every auxiliary atom below

    auto dense_d = std::vector<double>(nao * nao, 0.0);

    density.to_dense(dense_d.data());

    // NOTE: one of these per spin. An open shell has two, with a number of
    // orbitals of its own in each, and the exchange of a given auxiliary function
    // is their sum. A closed shell has one and reads exactly as it did.

    struct TDenseSpin
    {
        size_t              norbs;
        std::vector<double> dense_c;
        std::vector<double> transposed;
    };

    // NOTE: a positive omega is what makes this a range separated gradient. It
    // carries a second set of B vectors, a second set of fitted densities and a
    // second derivative tensor; what it does not carry is a second Coulomb term,
    // the attenuated operator entering the Fock matrix through the exchange alone.

    const auto range_separated = (omega > 0.0);

    errors::assertMsgCritical(range_separated || spins_erf.empty(),
                              std::string("SimdRIJKGradientDriver: Attenuated spins were given without a range "
                                          "separation parameter"));

    auto dense_spins = std::vector<TDenseSpin>();

    for (const auto &spin : spins)
    {
        const auto norbs = spin.coefficients->number_of_columns();

        auto dense_c = std::vector<double>(nao * norbs, 0.0);

        spin.coefficients->to_dense(dense_c.data());

        auto transposed = std::vector<double>(norbs * nao, 0.0);

        for (size_t mu = 0; mu < nao; mu++)
        {
            for (size_t ii = 0; ii < norbs; ii++) transposed[ii * nao + mu] = dense_c[mu * norbs + ii];
        }

        dense_spins.push_back({norbs, std::move(dense_c), std::move(transposed)});
    }

    // NOTE: the attenuated operator's spins carry the same orbitals -- a spin
    // occupies what it occupies whatever operator is being fitted -- so these are
    // the same dense arrays. They are built again rather than shared because the
    // list may be empty, and a reference into an empty list is worse than a copy.

    auto dense_spins_erf = std::vector<TDenseSpin>();

    for (const auto &spin : spins_erf)
    {
        const auto norbs = spin.coefficients->number_of_columns();

        auto dense_c = std::vector<double>(nao * norbs, 0.0);

        spin.coefficients->to_dense(dense_c.data());

        auto transposed = std::vector<double>(norbs * nao, 0.0);

        for (size_t mu = 0; mu < nao; mu++)
        {
            for (size_t ii = 0; ii < norbs; ii++) transposed[ii * nao + mu] = dense_c[mu * norbs + ii];
        }

        dense_spins_erf.push_back({norbs, std::move(dense_c), std::move(transposed)});
    }

    const auto naux = aux_basis.dimensions_of_basis();

    const auto aux_sets = aux_basis.basis_sets_indices();

    // every atom of the auxiliary basis, or the share the caller holds

    auto atoms = aux_atoms;

    if (atoms.empty())
    {
        atoms.resize(molecule.number_of_atoms());

        std::iota(atoms.begin(), atoms.end(), 0);
    }

    const auto eri_drv = CSimdThreeCenterElectronRepulsionDriver();

    const auto grad_drv = CSimdThreeCenterElectronRepulsionGradientDriver(_block_size);

    const auto grad_rs_drv = CSimdThreeCenterElectronRepulsionGradientRsDriver(_block_size);

    for (const auto iaux : atoms)
    {
        // NOTE: the pattern of this atom alone, described with the threshold the
        // calculation formed the B vectors with, so the derivative indexes the
        // same atom pairs they were formed over.

        const auto pattern = eri_drv.make_pattern(molecule, basis, aux_basis, _threshold, {iaux});

        // NOTE: the two operators come out of one call, so the attenuated term
        // costs a second contraction and not a second pass over the integrals.

        auto derivatives = std::vector<CSparseTensor>();

        if (range_separated)
        {
            auto pair = grad_rs_drv.compute(pattern, molecule, basis, aux_basis, omega);

            derivatives.push_back(std::move(pair.first));

            derivatives.push_back(std::move(pair.second));
        }
        else
        {
            derivatives.push_back(grad_drv.compute(pattern, molecule, basis, aux_basis));
        }

        // NOTE: the exchange part of Gamma, back transformed into the atomic
        // orbitals for every auxiliary function of this atom and held while its
        // derivative integrals are contracted. It was formed inside the loop over
        // the atom pairs before, which is the same matrix built again for every
        // pair of every block: two products of the basis functions squared by the
        // orbitals, where the whole term is meant to cost two of them per
        // auxiliary function.

        auto build_exchanges = [&](const std::vector<TExchangeSpin> &which,
                                   const std::vector<TDenseSpin>    &dense_which,
                                   const double                      factor,
                                   std::vector<size_t>              &slot_of,
                                   std::vector<std::vector<double>> &exchanges) {
            slot_of.assign(naux, naux);

            exchanges.clear();

            if (factor == 0.0) return;

            const auto &aux_set = aux_indices[static_cast<size_t>(aux_sets[static_cast<size_t>(iaux)])];

            for (const auto [lq, kq] : aux_set)
            {
                const auto nq = static_cast<size_t>(tensor::number_of_spherical_components(std::array<int, 1>{lq}));

                for (size_t mq = 0; mq < nq; mq++)
                {
                    const auto q = aux_starts[static_cast<size_t>(iaux) * aux_nmoms + lq] + kq + mq * aux_strides[lq];

                    auto matrix = std::vector<double>(nao * nao, 0.0);

                    // NOTE: the spins add into one matrix. For a closed shell the
                    // list is one long and this is the term it always was.

                    for (size_t is = 0; is < which.size(); is++)
                    {
                        const auto &dense = dense_which[is];

                        const auto norbs = dense.norbs;

                        if (norbs == 0) continue;

                        auto dense_q = std::vector<double>(norbs * norbs, 0.0);

                        (*which[is].orbital_densities)[q].to_dense(dense_q.data());

                        // C_o d(q), then that against the transposed orbitals

                        auto half = std::vector<double>(nao * norbs, 0.0);

                        _multiply(nao, norbs, norbs, dense.dense_c.data(), dense_q.data(), half.data());

                        // NOTE: added and not assigned. The spins accumulate into
                        // one matrix, and the overwriting product kept only the
                        // last of them.
                        _add_multiply(nao, nao, norbs, half.data(), dense.transposed.data(), matrix.data());
                    }

                    for (auto &value : matrix) value *= factor;

                    slot_of[q] = exchanges.size();

                    exchanges.push_back(std::move(matrix));
                }
            }
        };

        auto slot_of = std::vector<size_t>();

        auto exchanges = std::vector<std::vector<double>>();

        auto slot_of_erf = std::vector<size_t>();

        auto exchanges_erf = std::vector<std::vector<double>>();

        build_exchanges(spins, dense_spins, exchange_factor, slot_of, exchanges);

        if (range_separated)
        {
            build_exchanges(spins_erf, dense_spins_erf, erf_exchange_factor, slot_of_erf, exchanges_erf);
        }

        // NOTE: one operator at a time, and the same body for both of them. What
        // differs between the two calls is the derivative tensor, the weighted
        // densities it is contracted with, and whether there is a Coulomb term:
        // the attenuated operator has none.

        auto contract = [&](const CSparseTensor                    &derivative,
                            const double                            coulomb_factor,
                            const double                            exchange_factor,
                            const std::vector<std::vector<double>> &exchanges,
                            const std::vector<size_t>              &slot_of) {
        const auto nblocks = static_cast<size_t>(pattern.number_of_blocks());

        for (size_t iblk = 0; iblk < nblocks; iblk++)
        {
            const auto &block = pattern.block(iblk);

            const auto npairs = block.number_of_pairs();

            const auto natoms = block.number_of_c_atoms();

            if ((npairs == 0) || (natoms == 0)) continue;

            const auto &a_atoms = block.a_atoms();

            const auto &b_atoms = block.b_atoms();

            const auto &c_atoms = block.c_atoms();

            const auto &a_index = indices[static_cast<size_t>(block.a_index())];

            const auto &b_index = indices[static_cast<size_t>(block.b_index())];

            const auto &c_index = aux_indices[static_cast<size_t>(block.c_index())];

            for (size_t i = 0; i < a_index.size(); i++)
            {
                for (size_t j = 0; j < b_index.size(); j++)
                {
                    for (size_t k = 0; k < c_index.size(); k++)
                    {
                        const auto [la, ia] = a_index[i];

                        const auto [lb, jb] = b_index[j];

                        const auto [lc, kc] = c_index[k];

                        // NOTE: the atom pairs this combination reaches, which is
                        // what the kernel wrote and is not the atom pairs of the
                        // block: a combination which the screening cut short holds
                        // fewer, and they are the leading ones. Sizing the run by
                        // the block instead read past the end of its values.

                        const auto reached = block.number_of_pairs(la, ia, lb, jb, lc, kc);

                        if (reached == 0) continue;

                        const auto na = static_cast<size_t>(tensor::number_of_spherical_components(std::array<int, 1>{la}));

                        const auto nb = static_cast<size_t>(tensor::number_of_spherical_components(std::array<int, 1>{lb}));

                        const auto nc = static_cast<size_t>(tensor::number_of_spherical_components(std::array<int, 1>{lc}));

                        const auto ncomps = na * nb * nc;

                        const auto *values = derivative.values(iblk, la, ia, lb, jb, lc, kc);

                        const auto run = natoms * reached;

                        // NOTE: the components of a combination are the slowest
                        // index, the three of the first center then the three of
                        // the second, and the angular components sit inside them.

                        for (size_t ma = 0; ma < na; ma++)
                        {
                            for (size_t mb = 0; mb < nb; mb++)
                            {
                                for (size_t mc = 0; mc < nc; mc++)
                                {
                                    const auto angular = (ma * nb + mb) * nc + mc;

                                    for (size_t n = 0; n < natoms; n++)
                                    {
                                        const auto catom = static_cast<size_t>(c_atoms[n]);

                                        const auto q = aux_starts[catom * aux_nmoms + lc] + kc + mc * aux_strides[lc];

                                        // Gamma of this auxiliary function: the
                                        // Coulomb part is a scalar times the
                                        // density and is never stored, the
                                        // exchange part is the back transform of
                                        // the fitted density of the orbitals.

                                        const auto cp = coulomb_factor * fitting[q];

                                        const double *exchange =
                                            (exchange_factor != 0.0) ? exchanges[slot_of[q]].data() : nullptr;

                                        for (size_t p = 0; p < reached; p++)
                                        {
                                            const auto aatom = static_cast<size_t>(a_atoms[p]);

                                            const auto batom = static_cast<size_t>(b_atoms[p]);

                                            const auto mu = starts[aatom * nmoms + la] + ia + ma * strides[la];

                                            const auto nu = starts[batom * nmoms + lb] + jb + mb * strides[lb];

                                            auto gamma = cp * dense_d[mu * nao + nu];

                                            if (exchange != nullptr) gamma += exchange[mu * nao + nu];

                                            // NOTE: the pattern holds each unordered
                                            // pair of atoms once, so a pair of two
                                            // different atoms stands for the term of
                                            // the sum with its two orbitals the other
                                            // way round as well. A pair of one atom
                                            // holds both orders already and is not
                                            // counted twice.

                                            if (aatom != batom) gamma *= 2.0;

                                            for (size_t c = 0; c < 3; c++)
                                            {
                                                const auto at = ((c * ncomps) + angular) * run + n * reached + p;

                                                const auto bt = (((3 + c) * ncomps) + angular) * run + n * reached + p;

                                                const auto ta = gamma * values[at];

                                                const auto tb = gamma * values[bt];

                                                if (wanted[aatom]) gradient.data()[gradient.index(aatom, c)] += ta;

                                                if (wanted[batom]) gradient.data()[gradient.index(batom, c)] += tb;

                                                // NOTE: the three derivatives of
                                                // an integral sum to zero, so the
                                                // auxiliary center takes what the
                                                // other two leave.

                                                if (wanted[catom]) gradient.data()[gradient.index(catom, c)] -= ta + tb;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        };

        contract(derivatives[0], coulomb_factor, exchange_factor, exchanges, slot_of);

        if (range_separated)
        {
            // NOTE: a Coulomb factor of zero, and the fitting of the attenuated
            // operator is empty for the same reason: it has no Coulomb term. Passing
            // anything else here would read a vector which was deliberately left
            // unfilled.

            contract(derivatives[1], 0.0, erf_exchange_factor, exchanges_erf, slot_of_erf);
        }
    }
}
