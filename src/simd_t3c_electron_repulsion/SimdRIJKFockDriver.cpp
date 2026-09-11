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


#include "SimdRIJKFockDriver.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>

#include "ErrorHandler.hpp"
#include "PackedLinearAlgebra.hpp"

#ifdef VLX_USE_MATHLIB
#include "MathLibrary.hpp"
#else
#include "Eigen/Dense"
#endif
#include "ScreeningFunc.hpp"
#include "SimdThreeCenterElectronRepulsionDriver.hpp"
#include "SimdT3CDistributor.hpp"
#include "SimdTwoCenterElectronRepulsionDriver.hpp"
#include "TripleSparsityPattern.hpp"

auto
CSimdRIJKFockDriver::required_memory(const CMolecule       &molecule,
                                     const CMolecularBasis &basis,
                                     const CMolecularBasis &aux_basis,
                                     const double           threshold) const -> size_t
{
    // NOTE: the pattern of the B vectors is the pattern of the three-center
    // integrals, as the metric is dense and the transformation of the auxiliary
    // side keeps every atom which survives.

    const auto pattern = CSimdThreeCenterElectronRepulsionDriver().make_pattern(molecule, basis, aux_basis, threshold);

    size_t nvalues = 0;

    for (size_t i = 0; i < static_cast<size_t>(pattern.number_of_blocks()); i++)
    {
        nvalues += pattern.block(i).number_of_elements();
    }

    return nvalues * sizeof(double);
}

auto
CSimdRIJKFockDriver::prepare(const CMolecule       &molecule,
                             const CMolecularBasis &basis,
                             const CMolecularBasis &aux_basis,
                             const double           threshold,
                             const size_t           memory_budget,
                             const double           metric_threshold,
                             const bool             use_inverse_square_root,
                             const rimode           mode) -> void
{
    const auto memory = required_memory(molecule, basis, aux_basis, threshold);

    // NOTE: the memory is answered from the sparsity pattern, before any integral
    // is computed, so a molecule whose B vectors do not fit is put on the direct
    // way at once rather than after the work of forming them.

    _mode = (mode == rimode::automatic) ? ((memory > memory_budget) ? rimode::direct : rimode::in_memory) : mode;

    _molecule = molecule;

    _basis = basis;

    _aux_basis = aux_basis;

    const auto two_center = CSimdTwoCenterElectronRepulsionDriver().compute(molecule, aux_basis);

    // NOTE: both forms of the metric close the resolution of the identity, and the
    // Cholesky factor costs an order of magnitude less, so it is tried first. A
    // fitting basis which is close to linearly dependent has none, and the square
    // root is inverted in its place, dropping the directions which carry nothing.

    // NOTE: the direct way solves with the factor rather than multiplying by its
    // inverse, so it is the factor which is kept. The inverted square root is
    // taken for a metric which has no factor, in either way, as the B vectors
    // formed with it close the same sum.

    if (_mode == rimode::direct)
    {
        _pattern = CSimdThreeCenterElectronRepulsionDriver().make_pattern(molecule, basis, aux_basis, threshold);

        if (!use_inverse_square_root)
        {
            try
            {
                _factor = packlin::cholesky_factor(two_center);

                _bq_vectors = CSparseTensor();

                _w_vectors.clear();

                _prepared = true;

                return;
            }
            catch (const std::runtime_error &)
            {
                errors::msg(std::string("RIJKFockDriver: The metric of the fitting basis has no Cholesky factor. The "
                                        "direct way needs one, so the B vectors are formed with the inverted square "
                                        "root of the metric instead and the calculation is held in memory."),
                            "Warning");

                _mode = rimode::in_memory;
            }
        }
        else
        {
            errors::msg(std::string("RIJKFockDriver: The direct way solves with the Cholesky factor of the metric, "
                                    "which the inverted square root is not, so the calculation is held in memory."),
                        "Warning");

            _mode = rimode::in_memory;
        }
    }

    if (use_inverse_square_root)
    {
        _metric = packlin::inverse_square_root(two_center, metric_threshold);
    }
    else
    {
        try
        {
            _metric = packlin::cholesky_inverse(two_center);
        }
        catch (const std::runtime_error &)
        {
            errors::msg(std::string("RIJKFockDriver: The metric of the fitting basis has no Cholesky factor, so its "
                                    "square root is inverted instead. This is a nearly linearly dependent fitting "
                                    "basis."),
                        "Warning");

            _metric = packlin::inverse_square_root(two_center, metric_threshold);
        }
    }

    _bq_vectors = _drv.compute_bq_vectors(molecule, basis, aux_basis, _metric, threshold);

    _w_vectors.clear();

    _prepared = true;
}

auto
CSimdRIJKFockDriver::compute(const CPackedMatrix &density,
                             const CPackedMatrix &coefficients,
                             const double         exchange_scaling_factor) -> CPackedMatrix
{
    errors::assertMsgCritical(_prepared, std::string("RIJKFockDriver: The driver has not been prepared"));

    if (_mode == rimode::direct) return _compute_direct(density, coefficients, exchange_scaling_factor);

    // NOTE: the density of a closed shell calculation is that of one spin, so the
    // Coulomb matrix enters twice and the exchange once, scaled by the fraction of
    // exact exchange the functional asks for.

    auto fock = _drv.compute_fock_matrix(_bq_vectors, _basis, _aux_basis, density);

    fock.scale(2.0);

    if (exchange_scaling_factor == 0.0) return fock;

    const auto nao = _basis.dimensions_of_basis();

    const auto naux = _aux_basis.dimensions_of_basis();

    const auto norbitals = coefficients.number_of_columns();

    errors::assertMsgCritical((coefficients.get_type() == mat_t::general) && (coefficients.number_of_rows() == nao),
                              std::string("RIJKFockDriver: The orbital coefficients do not match the molecular basis"));

    if (norbitals == 0) return fock;

    // NOTE: the W matrices of a range are formed into storage the driver keeps, so
    // that the ranges of a call and the calls of a calculation reuse it. The number
    // of orbitals is taken from the coefficients rather than asked for, and the
    // storage is formed again only when it changes.

    const auto nbatch = std::min(_w_batch, naux);

    if ((_w_vectors.size() != nbatch) || (_w_vectors.front().number_of_columns() != norbitals) ||
        (_w_vectors.front().number_of_rows() != nao))
    {
        _w_vectors.clear();

        _w_vectors.reserve(nbatch);

        for (size_t i = 0; i < nbatch; i++)
        {
            _w_vectors.emplace_back(nao, norbitals, mat_t::general);
        }
    }

    for (size_t first = 0; first < naux; first += nbatch)
    {
        const auto last = std::min(first + nbatch, naux);

        const auto count = last - first;

        // NOTE: the last range is shorter than the others, and the storage is
        // handed to the transformation as the range it is asked to fill.

        if (count == nbatch)
        {
            _drv.compute_w_vectors(_bq_vectors, _basis, _aux_basis, coefficients, first, last, _w_vectors);

            _drv.compute_exchange_matrix(_w_vectors, fock, -exchange_scaling_factor);
        }
        else
        {
            auto tail = std::vector<CPackedMatrix>(_w_vectors.begin(), _w_vectors.begin() + static_cast<long>(count));

            _drv.compute_w_vectors(_bq_vectors, _basis, _aux_basis, coefficients, first, last, tail);

            _drv.compute_exchange_matrix(tail, fock, -exchange_scaling_factor);
        }
    }

    return fock;
}

auto
CSimdRIJKFockDriver::is_prepared() const -> bool
{
    return _prepared;
}

auto
CSimdRIJKFockDriver::get_mode() const -> rimode
{
    return _mode;
}

auto
CSimdRIJKFockDriver::get_bq_vectors() const -> const CSparseTensor &
{
    return _bq_vectors;
}

auto
CSimdRIJKFockDriver::get_metric() const -> const CPackedMatrix &
{
    return _metric;
}

auto
CSimdRIJKFockDriver::_compute_direct(const CPackedMatrix &density,
                                     const CPackedMatrix &coefficients,
                                     const double         exchange_scaling_factor) -> CPackedMatrix
{
    const auto nao = _basis.dimensions_of_basis();

    const auto naux = _aux_basis.dimensions_of_basis();

    const auto norbitals = coefficients.number_of_columns();

    errors::assertMsgCritical((coefficients.get_type() == mat_t::general) && (coefficients.number_of_rows() == nao),
                              std::string("RIJKFockDriver: The orbital coefficients do not match the molecular basis"));

    // NOTE: the half transformed integrals are the auxiliary basis by the basis
    // functions by the orbitals of a batch, and are held twice: once as the
    // matrices the transformation writes and once as the square the triangular
    // solve reads. The batch of orbitals is chosen so that the two fit.

    const auto per_orbital = 2 * naux * nao * sizeof(double);

    const auto nbatch = std::max(size_t{1}, std::min(norbitals, _direct_budget / std::max(per_orbital, size_t{1})));

    auto fock = CPackedMatrix(nao, nao, mat_t::symmetric);

    fock.zero();

    if (norbitals == 0) return fock;

    const auto *cvalues = coefficients.data();

    std::vector<double> gamma(naux, 0.0);

    CSimdThreeCenterElectronRepulsionDriver eri_drv;

    // the blocks of the pattern are formed in batches, so that the integrals of a
    // batch are held and dropped rather than all of them at once

    const auto &blocks = _pattern.blocks();

    std::vector<size_t> starts;

    {
        size_t memory = 0, first = 0;

        for (size_t i = 0; i < blocks.size(); i++)
        {
            const auto values = blocks[i].number_of_elements() * sizeof(double);

            if ((i > first) && (memory + values > _direct_budget))
            {
                starts.push_back(first);

                first = i;

                memory = 0;
            }

            memory += values;
        }

        starts.push_back(first);
    }

    starts.push_back(blocks.size());

    // the first pass: one sweep of the integrals for every batch of orbitals

    for (size_t ofirst = 0; ofirst < norbitals; ofirst += nbatch)
    {
        const auto olast = std::min(ofirst + nbatch, norbitals);

        const auto ncols = olast - ofirst;

        auto batch = CPackedMatrix(nao, ncols, mat_t::general);

        auto *bvalues = batch.data();

        for (size_t irow = 0; irow < nao; irow++)
        {
            std::copy(cvalues + irow * norbitals + ofirst, cvalues + irow * norbitals + olast,
                      bvalues + irow * ncols);
        }

        std::vector<CPackedMatrix> half;

        half.reserve(naux);

        for (size_t q = 0; q < naux; q++)
        {
            half.emplace_back(nao, ncols, mat_t::general);

            half.back().zero();
        }

        // NOTE: the half transformed integrals of one batch of orbitals are the
        // sum over every block of atom pairs, so the blocks are swept and added
        // into rather than each one setting them.

        for (size_t k = 0; k + 1 < starts.size(); k++)
        {
            auto part = std::vector<CAtomBasisTripleSparsity>(blocks.begin() + static_cast<long>(starts[k]),
                                                              blocks.begin() + static_cast<long>(starts[k + 1]));

            if (part.empty()) continue;

            const auto pattern = CTripleSparsityPattern(std::move(part), mat_t::symmetric, _pattern.get_threshold());

            auto integrals = CSparseTensor(pattern);

            integrals.allocate();

            auto distributor = CSimdT3CDistributor<CSparseTensor>(&integrals);

            eri_drv.compute(pattern, _molecule, _basis, _aux_basis, distributor);

            _drv.compute_w_vectors(integrals, _basis, _aux_basis, batch, 0, naux, half, true);
        }

        // the Coulomb vector, which is the half transformed integrals closed with
        // the orbitals they were transformed by

        for (size_t q = 0; q < naux; q++)
        {
            const auto *values = half[q].data();

            double sum = 0.0;

            for (size_t irow = 0; irow < nao; irow++)
            {
                const auto *crow = bvalues + irow * ncols;

                const auto *vrow = values + irow * ncols;

#pragma omp simd reduction(+ : sum)
                for (size_t j = 0; j < ncols; j++)
                {
                    sum += vrow[j] * crow[j];
                }
            }

            gamma[q] += sum;
        }

        // the exchange: solve the factor against the half transformed integrals,
        // which gives the B vectors of this batch of orbitals, and add their
        // square into the matrix

        const auto width = nao * ncols;

        std::vector<double> stacked(naux * width);

        for (size_t q = 0; q < naux; q++)
        {
            std::copy(half[q].data(), half[q].data() + width, stacked.data() + q * width);
        }

        _solve_factor(stacked.data(), naux, width, false);

        for (size_t q = 0; q < naux; q++)
        {
            std::copy(stacked.data() + q * width, stacked.data() + (q + 1) * width, half[q].data());
        }

        _drv.compute_exchange_matrix(half, fock, -exchange_scaling_factor);
    }

    // the coefficients of the fitting, from the factor and its transpose

    _solve_factor(gamma.data(), naux, 1, false);

    _solve_factor(gamma.data(), naux, 1, true);

    // the second pass: the Coulomb matrix from the integrals and those coefficients

    for (size_t k = 0; k + 1 < starts.size(); k++)
    {
        auto part = std::vector<CAtomBasisTripleSparsity>(blocks.begin() + static_cast<long>(starts[k]),
                                                          blocks.begin() + static_cast<long>(starts[k + 1]));

        if (part.empty()) continue;

        const auto pattern = CTripleSparsityPattern(std::move(part), mat_t::symmetric, _pattern.get_threshold());

        auto integrals = CSparseTensor(pattern);

        integrals.allocate();

        auto distributor = CSimdT3CDistributor<CSparseTensor>(&integrals);

        eri_drv.compute(pattern, _molecule, _basis, _aux_basis, distributor);

        // NOTE: the blocks of atom pairs of one batch write elements of the matrix
        // no other batch writes, so the parts are added and nothing is counted
        // twice. The Coulomb matrix enters twice, as the density is of one spin.

        auto part_fock = _drv.compute_fock_matrix(integrals, _basis, _aux_basis, gamma);

        auto       *values = fock.data();

        const auto *added = part_fock.data();

        const auto nvalues = static_cast<int>(fock.number_of_elements());

#pragma omp parallel for schedule(static) if (nvalues > 1)
        for (int i = 0; i < nvalues; i++)
        {
            values[static_cast<size_t>(i)] += 2.0 * added[static_cast<size_t>(i)];
        }
    }

    return fock;
}

auto
CSimdRIJKFockDriver::_solve_factor(double *values, const size_t nrows, const size_t ncols, const bool transposed) const
    -> void
{
    // NOTE: the factor is held in the packed format and the library solves against
    // a square, so it is expanded once here. It is the square of the auxiliary
    // basis, which is small beside the half transformed integrals it solves.

    auto dense = std::vector<double>(nrows * nrows, 0.0);

    _factor.to_dense(dense.data());

#ifdef VLX_USE_MATHLIB

    // NOTE: the array of the factor is row major and the library is column major,
    // so the column major matrix of it is the transpose of the lower triangular
    // factor, which is upper triangular. Solving L X = B with the row major arrays
    // is therefore the transposed upper triangular solve of the library, and the
    // right hand sides are transposed in the same way, which swaps the side.

    const char side = 'R';

    const char uplo = 'U';

    // NOTE: solving L X = B with row major arrays is X L transposed = B
    // transposed, which the library takes as a solve from the right against the
    // upper triangular matrix its column major reading of the factor gives. That
    // is the untransposed call; solving against the transpose of the factor is the
    // transposed one.

    const char trans = transposed ? 'T' : 'N';

    const char diag = 'N';

    auto m_arg = static_cast<lapack_int_t>(ncols);

    auto n_arg = static_cast<lapack_int_t>(nrows);

    auto lda = static_cast<lapack_int_t>(nrows);

    auto ldb = static_cast<lapack_int_t>(ncols);

    const double one = 1.0;

    dtrsm_(&side, &uplo, &trans, &diag, &m_arg, &n_arg, &one, dense.data(), &lda, values, &ldb);

#else

    using RowMajorMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    const auto n = static_cast<Eigen::Index>(nrows);

    const auto m = static_cast<Eigen::Index>(ncols);

    Eigen::Map<const RowMajorMatrix> lmap(dense.data(), n, n);

    Eigen::Map<RowMajorMatrix> bmap(values, n, m);

    if (transposed)
    {
        lmap.transpose().template triangularView<Eigen::Upper>().solveInPlace(bmap);
    }
    else
    {
        lmap.template triangularView<Eigen::Lower>().solveInPlace(bmap);
    }

#endif /* VLX_USE_MATHLIB */
}
