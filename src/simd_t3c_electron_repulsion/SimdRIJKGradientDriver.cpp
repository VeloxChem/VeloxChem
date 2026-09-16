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

#include "DenseIndexFunc.hpp"
#include "ErrorHandler.hpp"
#include "SimdThreeCenterElectronRepulsionDriver.hpp"
#include "SimdThreeCenterElectronRepulsionGradientDriver.hpp"
#include "SimdTwoCenterElectronRepulsionGradientDriver.hpp"
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

/// @brief Applies the transpose of the inverted Cholesky factor to a set of
/// columns held one after another.
/// @note The factor is held as the inverted lower triangular matrix, so its
/// transpose is upper triangular and no second matrix is stored for it: the
/// transposed multiply is the same array read the other way.
/// @brief The same, over a factor which the caller has already expanded.
///
/// @note The expansion is the square of the auxiliary basis and is the same array
/// every time. Forming it inside this routine costs nothing when the routine is
/// called once, and everything when it is called once per element of a matrix:
/// ninety-two megabytes written eight thousand times over on a molecule of this
/// size, which is memory traffic of a wholly different order to the arithmetic it
/// serves. The caller which loops expands it once and hands it in.
static auto
_multiply_transposed(const double *factor, const size_t naux, double *values, const size_t ncols) -> void
{
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

static auto
_multiply_transposed(const CPackedMatrix &metric, double *values, const size_t ncols) -> void
{
    const auto naux = metric.number_of_rows();

    auto factor = std::vector<double>(naux * naux, 0.0);

    metric.to_dense(factor.data());

    _multiply_transposed(factor.data(), naux, values, ncols);
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

    // NOTE: the transpose of the factor, written out once. The element of the
    // result at p reads the rows at and below p, which is the transpose of a
    // lower triangular matrix read as it stands, and a general multiply wants it
    // that way round. The expansion it is built from is released before the
    // arrays which are the size of the whole phase are taken.

    auto transposed = std::vector<double>(naux * naux, 0.0);

    {
        auto factor = std::vector<double>(naux * naux, 0.0);

        metric.to_dense(factor.data());

        for (size_t p = 0; p < naux; p++)
        {
            for (size_t q = p; q < naux; q++)
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

    auto gathered = std::vector<double>(naux * nelements, 0.0);

    const auto nrange = static_cast<int>(naux);

#pragma omp parallel for schedule(static)
    for (int at = 0; at < nrange; at++)
    {
        const auto q = static_cast<size_t>(at);

        std::copy_n(matrices[q].data(), nelements, gathered.data() + q * nelements);
    }

    auto result = std::vector<double>(naux * nelements, 0.0);

    _multiply(naux, nelements, naux, transposed.data(), gathered.data(), result.data());

#pragma omp parallel for schedule(static)
    for (int at = 0; at < nrange; at++)
    {
        const auto q = static_cast<size_t>(at);

        std::copy_n(result.data() + q * nelements, nelements, matrices[q].data());
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

    auto wanted = std::vector<bool>(natoms, false);

    for (const auto iatom : atoms) wanted[static_cast<size_t>(iatom)] = true;

    // the fitted densities, which read no integrals

    const auto fitted = fitted_densities(bq_vectors, basis, aux_basis, metric, density, coefficients,
                                         exchange_scaling_factor);

    // the three-center term, plus in the gradient

    _compute_three_center(gradient, molecule, basis, aux_basis, fitted, density, coefficients,
                          exchange_scaling_factor, wanted, aux_atoms);

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
CSimdRIJKGradientDriver::_compute_three_center(CPackedMatrix           &gradient,
                                               const CMolecule         &molecule,
                                               const CMolecularBasis   &basis,
                                               const CMolecularBasis   &aux_basis,
                                               const TFittedDensities  &fitted,
                                               const CPackedMatrix     &density,
                                               const CPackedMatrix     &coefficients,
                                               const double             exchange_scaling_factor,
                                               const std::vector<bool> &wanted,
                                               const std::vector<int>  &aux_atoms) const -> void
{
    const auto nao = coefficients.number_of_rows();

    const auto norbs = coefficients.number_of_columns();

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

    auto dense_c = std::vector<double>(nao * norbs, 0.0);

    coefficients.to_dense(dense_c.data());

    auto transposed = std::vector<double>(norbs * nao, 0.0);

    for (size_t mu = 0; mu < nao; mu++)
    {
        for (size_t ii = 0; ii < norbs; ii++) transposed[ii * nao + mu] = dense_c[mu * norbs + ii];
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

    for (const auto iaux : atoms)
    {
        // NOTE: the pattern of this atom alone, described with the threshold the
        // calculation formed the B vectors with, so the derivative indexes the
        // same atom pairs they were formed over.

        const auto pattern = eri_drv.make_pattern(molecule, basis, aux_basis, _threshold, {iaux});

        const auto derivative = grad_drv.compute(pattern, molecule, basis, aux_basis);

        // NOTE: the exchange part of Gamma, back transformed into the atomic
        // orbitals for every auxiliary function of this atom and held while its
        // derivative integrals are contracted. It was formed inside the loop over
        // the atom pairs before, which is the same matrix built again for every
        // pair of every block: two products of the basis functions squared by the
        // orbitals, where the whole term is meant to cost two of them per
        // auxiliary function.

        auto slot_of = std::vector<size_t>(naux, naux);

        auto exchanges = std::vector<std::vector<double>>();

        if (exchange_scaling_factor != 0.0)
        {
            const auto &aux_set = aux_indices[static_cast<size_t>(aux_sets[static_cast<size_t>(iaux)])];

            auto half = std::vector<double>(nao * norbs, 0.0);

            auto dense_q = std::vector<double>(norbs * norbs, 0.0);

            for (const auto [lq, kq] : aux_set)
            {
                const auto nq = static_cast<size_t>(tensor::number_of_spherical_components(std::array<int, 1>{lq}));

                for (size_t mq = 0; mq < nq; mq++)
                {
                    const auto q = aux_starts[static_cast<size_t>(iaux) * aux_nmoms + lq] + kq + mq * aux_strides[lq];

                    fitted.orbital_densities[q].to_dense(dense_q.data());

                    // C_o d(q), then that against the transposed orbitals

                    _multiply(nao, norbs, norbs, dense_c.data(), dense_q.data(), half.data());

                    auto matrix = std::vector<double>(nao * nao, 0.0);

                    _multiply(nao, nao, norbs, half.data(), transposed.data(), matrix.data());

                    for (auto &value : matrix) value *= -2.0 * exchange_scaling_factor;

                    slot_of[q] = exchanges.size();

                    exchanges.push_back(std::move(matrix));
                }
            }
        }

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

                                        const auto cp = 4.0 * fitted.coefficients[q];

                                        const double *exchange =
                                            (exchange_scaling_factor != 0.0) ? exchanges[slot_of[q]].data() : nullptr;

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
    }
}
