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
#include <numeric>
#include <set>
#include <string>

#include "ErrorHandler.hpp"

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

/// @brief Applies the transpose of the inverted Cholesky factor to a set of
/// columns held one after another.
/// @note The factor is held as the inverted lower triangular matrix, so its
/// transpose is upper triangular and no second matrix is stored for it: the
/// transposed multiply is the same array read the other way.
static auto
_multiply_transposed(const CPackedMatrix &metric, double *values, const size_t ncols) -> void
{
    const auto naux = metric.number_of_rows();

    auto factor = std::vector<double>(naux * naux, 0.0);

    metric.to_dense(factor.data());

    auto column = std::vector<double>(naux, 0.0);

    for (size_t icol = 0; icol < ncols; icol++)
    {
        auto *target = values + icol * naux;

        std::fill(column.begin(), column.end(), 0.0);

        // NOTE: the transpose of a lower triangular matrix, so the element p of
        // the result reads the rows at and below p of the column it multiplies.

        for (size_t p = 0; p < naux; p++)
        {
            double sum = 0.0;

            for (size_t q = p; q < naux; q++)
            {
                sum += factor[q * naux + p] * target[q];
            }

            column[p] = sum;
        }

        std::copy(column.begin(), column.end(), target);
    }
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

    auto column = std::vector<double>(naux, 0.0);

    for (size_t at = 0; at < nelements; at++)
    {
        for (size_t q = 0; q < naux; q++)
        {
            column[q] = matrices[q].data()[at];
        }

        _multiply_transposed(metric, column.data(), 1);

        for (size_t q = 0; q < naux; q++)
        {
            matrices[q].data()[at] = column[q];
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
    // the factor and once as its transpose, which the inverted Cholesky factor
    // gives and the inverted square root does not: the two close the same sum for
    // the energy and are not the same matrix. Which of them it was handed is read
    // from the type -- a factor is lower triangular and a root is symmetric --
    // rather than from a flag beside it, the way the Fock driver reads it.

    errors::assertMsgCritical(
        metric.number_of_rows() > 0,
        std::string("SimdRIJKGradientDriver: The metric is empty. The direct mode forms no inverted "
                    "factor, and this driver needs the one the mode which holds the B vectors forms"));

    errors::assertMsgCritical(
        metric.get_type() != mat_t::symmetric,
        std::string("SimdRIJKGradientDriver: The metric is the inverted square root, which this driver "
                    "does not support. Set ri_metric_route to cholesky for a gradient"));

    errors::assertMsgCritical(
        metric.number_of_rows() == aux_basis.dimensions_of_basis(),
        std::string("SimdRIJKGradientDriver: The metric is not of the auxiliary basis"));
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

    const auto nao = coefficients.number_of_rows();

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

    // and the two index fitted density, which the derivative of the metric is
    // contracted against. The exchange part is the Gram product of the fitted
    // densities of the orbitals and the Coulomb part is the outer product of the
    // fitting coefficients.

    auto omega = CPackedMatrix(naux, naux, mat_t::symmetric);

    omega.zero();

    for (size_t p = 0; p < naux; p++)
    {
        for (size_t q = 0; q <= p; q++)
        {
            auto value = 2.0 * fitting[p] * fitting[q];

            if (exchange_scaling_factor != 0.0)
            {
                const auto *dp = orbital_densities[p].data();

                const auto *dq = orbital_densities[q].data();

                double sum = 0.0;

                // NOTE: the packed matrices hold the lower triangle, so the
                // elements off the diagonal stand for two of the sum over the
                // pairs of orbitals and are counted twice.

                for (size_t i = 0; i < norbs; i++)
                {
                    for (size_t j = 0; j < i; j++)
                    {
                        sum += 2.0 * dp[i * (i + 1) / 2 + j] * dq[i * (i + 1) / 2 + j];
                    }

                    sum += dp[i * (i + 1) / 2 + i] * dq[i * (i + 1) / 2 + i];
                }

                value -= exchange_scaling_factor * sum;
            }

            omega.data()[omega.index(p, q)] = value;
        }
    }

    return {std::move(fitting), std::move(orbital_densities), std::move(omega)};
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

    // TODO: the derivative integrals and their contraction. The kernels which
    // form them are not written yet -- the plain tree has the derivative of the
    // three-center integral with respect to the auxiliary center alone, and none
    // of them are vectorised -- so this driver is a skeleton whose arguments are
    // fixed and whose body is not.
    //
    // What it will do, once they are: form the fitting coefficients from the B
    // vectors and the density, then for each atom of `atoms` accumulate the
    // derivative of the three-center integrals against them for the Coulomb part
    // and against the W matrices of the occupied orbitals for the exchange part,
    // and the derivative of the metric against the coefficients on both sides.
    //
    // The unused arguments are named and not commented out so that the shape of
    // the call does not change when the body arrives.
    // Phase one, which is written: the fitted densities. Everything below it
    // waits on the derivative integrals.
    const auto fitted = fitted_densities(bq_vectors, basis, aux_basis, metric, density, coefficients,
                                         exchange_scaling_factor);

    (void)fitted;
    (void)aux_atoms;

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
