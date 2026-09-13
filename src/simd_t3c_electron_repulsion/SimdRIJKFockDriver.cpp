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
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <memory>
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
#include "OpenMPFunc.hpp"
#include "TripleSparsityPattern.hpp"

namespace {

/// @brief The clock the phases of a Fock build are timed on.
using prof_clock = std::chrono::steady_clock;

/// @brief Whether the phases of a Fock build are to be reported.
/// @note Answered once, from the environment, so that an ordinary calculation pays
/// nothing for the question.
inline auto
prof_wanted() -> bool
{
    static const bool wanted = (std::getenv("VLX_RIJK_PROFILE") != nullptr);

    return wanted;
}

/// @brief The seconds since a mark.
inline auto
prof_since(const prof_clock::time_point &mark) -> double
{
    return std::chrono::duration<double>(prof_clock::now() - mark).count();
}

/// @brief The times of the phases of one Fock build of the direct mode.
/// @note The direct mode repeats every phase on every iteration, so a phase which
/// does not widen with the threads bounds the whole calculation however many cores
/// are given to it. Run a calculation at one thread and again at many, with
/// VLX_RIJK_PROFILE set, and the phase whose time does not fall between the two is
/// the one worth working on. The whole is timed as well as the parts, so that what
/// the parts do not account for is visible rather than assumed.
struct CDirectProfile
{
    double allocate = 0.0;
    double integrals_a = 0.0;
    double transform = 0.0;
    double closure = 0.0;
    double copies = 0.0;
    double solve = 0.0;
    double exchange = 0.0;
    double integrals_b = 0.0;
    double coulomb = 0.0;
    double total = 0.0;

    /// @brief Writes the phases of this build, and their share of it.
    auto report() const -> void
    {
        static size_t builds = 0;

        builds++;

        const auto accounted = allocate + integrals_a + transform + closure + copies +
                               solve + exchange + integrals_b + coulomb;

        const char *names[] = {"allocate", "integrals a", "transform", "closure", "copies",
                               "solve",    "exchange",    "integrals b", "coulomb", "rest"};

        const double times[] = {allocate, integrals_a, transform, closure,  copies,
                                solve,    exchange,    integrals_b, coulomb, total - accounted};

        std::printf("RIJK build %zu on %d threads, %.3f s\n", builds, omp::get_number_of_threads(), total);

        for (size_t i = 0; i < 10; i++)
        {
            std::printf("RIJK   %-12s %9.3f s %6.1f %%\n", names[i], times[i],
                        (total > 0.0) ? 100.0 * times[i] / total : 0.0);
        }

        std::fflush(stdout);
    }
};

/// @brief The times of the phases of one Fock build of the mode which holds the B
/// vectors.
/// @note The phases are those of the direct mode which this one does not repeat.
/// It forms no integrals, so what is left is the Coulomb matrix, the transformation
/// which forms the W matrices of a range, and the exchange. The exchange writes a
/// line of its own from the driver beneath this one, and is counted here as well so
/// that what the phases leave over is the rest and nothing else.
struct CInMemoryProfile
{
    double coulomb = 0.0;
    double allocate = 0.0;
    double transform = 0.0;
    double exchange = 0.0;
    double total = 0.0;

    /// @brief Writes the phases of this build, and their share of it.
    auto report() const -> void
    {
        static size_t builds = 0;

        builds++;

        const auto accounted = coulomb + allocate + transform + exchange;

        const char *names[] = {"coulomb", "allocate", "transform", "exchange", "rest"};

        const double times[] = {coulomb, allocate, transform, exchange, total - accounted};

        std::printf("RIJK memory build %zu on %d threads, %.3f s\n", builds, omp::get_number_of_threads(), total);

        for (size_t i = 0; i < 5; i++)
        {
            std::printf("RIJK   %-12s %9.3f s %6.1f %%\n", names[i], times[i],
                        (total > 0.0) ? 100.0 * times[i] / total : 0.0);
        }

        std::fflush(stdout);
    }
};

/// @brief The times of the phases of the setup of a calculation.
/// @note The setup is paid once, where the phases of a build are paid on every
/// iteration, so the two reports are not to be added. It is here because the mode
/// which holds the B vectors forms them here and the direct mode does not, which is
/// the whole of what the two modes do differently outside a Fock build.
struct CPrepareProfile
{
    double two_center = 0.0;
    double pattern = 0.0;
    double metric = 0.0;
    double bq_vectors = 0.0;
    double total = 0.0;

    /// @brief Writes the phases of the setup, and their share of it.
    auto report(const char *mode) const -> void
    {
        const auto accounted = two_center + pattern + metric + bq_vectors;

        const char *names[] = {"two center", "pattern", "metric", "b vectors", "rest"};

        const double times[] = {two_center, pattern, metric, bq_vectors, total - accounted};

        std::printf("RIJK setup of the %s way on %d threads, %.3f s\n", mode, omp::get_number_of_threads(), total);

        for (size_t i = 0; i < 5; i++)
        {
            std::printf("RIJK   %-12s %9.3f s %6.1f %%\n", names[i], times[i],
                        (total > 0.0) ? 100.0 * times[i] / total : 0.0);
        }

        std::fflush(stdout);
    }
};

}  // namespace

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
    CPrepareProfile profile;

    const auto profile_start = prof_clock::now();

    const auto memory = required_memory(molecule, basis, aux_basis, threshold);

    _budget = memory_budget;

    // NOTE: the memory is answered from the sparsity pattern, before any integral
    // is computed, so a molecule whose B vectors do not fit is put on the direct
    // way at once rather than after the work of forming them.

    _mode = (mode == rimode::automatic) ? ((memory > memory_budget) ? rimode::direct : rimode::in_memory) : mode;

    _molecule = molecule;

    _basis = basis;

    _aux_basis = aux_basis;

    const auto mark_two_center = prof_clock::now();

    const auto two_center = CSimdTwoCenterElectronRepulsionDriver().compute(molecule, aux_basis);

    profile.two_center += prof_since(mark_two_center);

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
        const auto mark_pattern = prof_clock::now();

        const auto pattern = CSimdThreeCenterElectronRepulsionDriver().make_pattern(molecule, basis, aux_basis, threshold);

        _parts = _make_parts(molecule, basis, aux_basis, threshold, pattern);

        profile.pattern += prof_since(mark_pattern);

        if (!use_inverse_square_root)
        {
            try
            {
                const auto mark_metric = prof_clock::now();

                _factor = packlin::cholesky_factor(two_center);

                profile.metric += prof_since(mark_metric);

                _bq_vectors = CSparseTensor();

                _w_vectors.clear();

                _prepared = true;

                profile.total = prof_since(profile_start);

                if (prof_wanted()) profile.report("direct");

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

    const auto mark_metric = prof_clock::now();

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

    profile.metric += prof_since(mark_metric);

    const auto mark_bq_vectors = prof_clock::now();

    _bq_vectors = _drv.compute_bq_vectors(molecule, basis, aux_basis, _metric, threshold);

    profile.bq_vectors += prof_since(mark_bq_vectors);

    _w_vectors.clear();

    _prepared = true;

    profile.total = prof_since(profile_start);

    if (prof_wanted()) profile.report("memory");
}

auto
CSimdRIJKFockDriver::compute(const CPackedMatrix &density,
                             const CPackedMatrix &coefficients,
                             const double         exchange_scaling_factor) -> CPackedMatrix
{
    errors::assertMsgCritical(_prepared, std::string("RIJKFockDriver: The driver has not been prepared"));

    if (_mode == rimode::direct) return _compute_direct(density, coefficients, exchange_scaling_factor);

    CInMemoryProfile profile;

    const auto profile_start = prof_clock::now();

    // NOTE: the density of a closed shell calculation is that of one spin, so the
    // Coulomb matrix enters twice and the exchange once, scaled by the fraction of
    // exact exchange the functional asks for.

    const auto mark_coulomb = prof_clock::now();

    auto fock = _drv.compute_fock_matrix(_bq_vectors, _basis, _aux_basis, density);

    fock.scale(2.0);

    profile.coulomb += prof_since(mark_coulomb);

    if (exchange_scaling_factor == 0.0)
    {
        profile.total = prof_since(profile_start);

        if (prof_wanted()) profile.report();

        return fock;
    }

    const auto nao = _basis.dimensions_of_basis();

    const auto naux = _aux_basis.dimensions_of_basis();

    const auto norbitals = coefficients.number_of_columns();

    errors::assertMsgCritical((coefficients.get_type() == mat_t::general) && (coefficients.number_of_rows() == nao),
                              std::string("RIJKFockDriver: The orbital coefficients do not match the molecular basis"));

    if (norbitals == 0)
    {
        profile.total = prof_since(profile_start);

        if (prof_wanted()) profile.report();

        return fock;
    }

    // NOTE: the W matrices of a range are formed into storage the driver keeps, so
    // that the ranges of a call and the calls of a calculation reuse it. The number
    // of orbitals is taken from the coefficients rather than asked for, and the
    // storage is formed again only when it changes.

    // NOTE: the range is as long as the memory allows. The transformation divides
    // it over the threads, so it has to be at least as long as there are threads,
    // but what it costs beside its tasks -- the parallel region, and a square of
    // the basis allocated and zeroed for every thread -- is paid once for the call
    // however long the range is. Taking the longest range the memory allows makes
    // the fewest calls, and each function of the range holds a matrix of the basis
    // by the occupied orbitals.

    const auto per_function = nao * norbitals * sizeof(double);

    const auto by_memory = std::max(size_t{1}, _w_batch_memory / std::max(per_function, size_t{1}));

    const auto nbatch = std::min(naux, std::max(_w_batch, by_memory));

    const auto mark_allocate = prof_clock::now();

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

    profile.allocate += prof_since(mark_allocate);

    for (size_t first = 0; first < naux; first += nbatch)
    {
        const auto last = std::min(first + nbatch, naux);

        const auto count = last - first;

        // NOTE: the last range is shorter than the others, and the storage is
        // handed to the transformation as the range it is asked to fill.

        const auto mark_transform = prof_clock::now();

        if (count == nbatch)
        {
            _drv.compute_w_vectors(_bq_vectors, _basis, _aux_basis, coefficients, first, last, _w_vectors);

            profile.transform += prof_since(mark_transform);

            const auto mark_exchange = prof_clock::now();

            _drv.compute_exchange_matrix(_w_vectors, fock, -exchange_scaling_factor);

            profile.exchange += prof_since(mark_exchange);
        }
        else
        {
            auto tail = std::vector<CPackedMatrix>(_w_vectors.begin(), _w_vectors.begin() + static_cast<long>(count));

            _drv.compute_w_vectors(_bq_vectors, _basis, _aux_basis, coefficients, first, last, tail);

            profile.transform += prof_since(mark_transform);

            const auto mark_exchange = prof_clock::now();

            _drv.compute_exchange_matrix(tail, fock, -exchange_scaling_factor);

            profile.exchange += prof_since(mark_exchange);
        }
    }

    profile.total = prof_since(profile_start);

    if (prof_wanted()) profile.report();

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

    const auto nbatch = std::max(size_t{1}, std::min(norbitals, (_budget / 2) / std::max(per_orbital, size_t{1})));

    CDirectProfile profile;

    const auto profile_start = prof_clock::now();

    auto fock = CPackedMatrix(nao, nao, mat_t::symmetric);

    fock.zero();

    if (norbitals == 0) return fock;

    const auto *cvalues = coefficients.data();

    std::vector<double> gamma(naux, 0.0);

    CSimdThreeCenterElectronRepulsionDriver eri_drv;

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

        // NOTE: the matrix of one auxiliary function is smaller than the chunk the
        // packed matrix divides its zeroing by, so its constructor zeroes it on the
        // calling thread alone. Several thousand of them is several gigabytes of
        // memset, and of first touch besides, which is why the functions are
        // divided over the threads here and each thread allocates and zeroes its
        // own. The constructor is what zeroes them; nothing zeroes them twice.

        const auto mark_allocate = prof_clock::now();

        std::vector<CPackedMatrix> half(naux);

        const auto nmatrices = static_cast<int>(naux);

#pragma omp parallel for schedule(static) if (nmatrices > 1)
        for (int iq = 0; iq < nmatrices; iq++)
        {
            half[static_cast<size_t>(iq)] = CPackedMatrix(nao, ncols, mat_t::general);
        }

        profile.allocate += prof_since(mark_allocate);

        // NOTE: the half transformed integrals of one batch of orbitals are the
        // sum over every block of atom pairs, so the blocks are swept and added
        // into rather than each one setting them.

        for (const auto &pattern : _parts)
        {
            auto integrals = CSparseTensor(pattern);

            integrals.allocate();

            auto distributor = CSimdT3CDistributor<CSparseTensor>(&integrals);

            const auto mark_integrals = prof_clock::now();

            eri_drv.compute(pattern, _molecule, _basis, _aux_basis, distributor);

            profile.integrals_a += prof_since(mark_integrals);

            const auto mark_transform = prof_clock::now();

            _drv.compute_w_vectors(integrals, _basis, _aux_basis, batch, 0, naux, half, true);

            profile.transform += prof_since(mark_transform);
        }

        // the Coulomb vector, which is the half transformed integrals closed with
        // the orbitals they were transformed by

        // NOTE: an auxiliary function is closed against the orbitals it was
        // transformed by, and each one adds into a place of its own, so the
        // functions are divided over the threads without anything to reduce.

        const auto mark_closure = prof_clock::now();

        const auto nrange = static_cast<int>(naux);

#pragma omp parallel for schedule(static) if (nrange > 1)
        for (int iq = 0; iq < nrange; iq++)
        {
            const auto q = static_cast<size_t>(iq);

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

        profile.closure += prof_since(mark_closure);

        // the exchange: solve the factor against the half transformed integrals,
        // which gives the B vectors of this batch of orbitals, and add their
        // square into the matrix

        const auto width = nao * ncols;

        // NOTE: every element of the array is written by the copy below before it is
        // read, so it is left with the content of the allocation rather than zeroed
        // first. Zeroing it would be a sweep of several gigabytes for nothing.

        auto mark_copies = prof_clock::now();

        auto stacked = std::make_unique_for_overwrite<double[]>(naux * width);

        const auto ncopies = static_cast<int>(naux);

#pragma omp parallel for schedule(static) if (ncopies > 1)
        for (int iq = 0; iq < ncopies; iq++)
        {
            const auto q = static_cast<size_t>(iq);

            std::copy(half[q].data(), half[q].data() + width, stacked.get() + q * width);
        }

        profile.copies += prof_since(mark_copies);

        const auto mark_solve = prof_clock::now();

        _solve_factor(stacked.get(), naux, width, false);

        profile.solve += prof_since(mark_solve);

        mark_copies = prof_clock::now();

#pragma omp parallel for schedule(static) if (ncopies > 1)
        for (int iq = 0; iq < ncopies; iq++)
        {
            const auto q = static_cast<size_t>(iq);

            std::copy(stacked.get() + q * width, stacked.get() + (q + 1) * width, half[q].data());
        }

        profile.copies += prof_since(mark_copies);

        const auto mark_exchange = prof_clock::now();

        _drv.compute_exchange_matrix(half, fock, -exchange_scaling_factor);

        profile.exchange += prof_since(mark_exchange);
    }

    // the coefficients of the fitting, from the factor and its transpose

    const auto mark_gamma = prof_clock::now();

    _solve_factor(gamma.data(), naux, 1, false);

    _solve_factor(gamma.data(), naux, 1, true);

    profile.solve += prof_since(mark_gamma);

    // the second pass: the Coulomb matrix from the integrals and those coefficients

    for (const auto &pattern : _parts)
    {
        auto integrals = CSparseTensor(pattern);

        integrals.allocate();

        auto distributor = CSimdT3CDistributor<CSparseTensor>(&integrals);

        const auto mark_integrals = prof_clock::now();

        eri_drv.compute(pattern, _molecule, _basis, _aux_basis, distributor);

        profile.integrals_b += prof_since(mark_integrals);

        const auto mark_coulomb = prof_clock::now();

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

        profile.coulomb += prof_since(mark_coulomb);
    }

    profile.total = prof_since(profile_start);

    if (prof_wanted()) profile.report();

    return fock;
}

auto
CSimdRIJKFockDriver::_make_parts(const CMolecule              &molecule,
                                 const CMolecularBasis        &basis,
                                 const CMolecularBasis        &aux_basis,
                                 const double                  threshold,
                                 const CTripleSparsityPattern &pattern) const -> std::vector<CTripleSparsityPattern>
{
    const auto natoms = static_cast<size_t>(molecule.number_of_atoms());

    // NOTE: a block holds as many values for one of its atoms on the auxiliary side
    // as for any other, so its memory divides evenly over them and the memory of an
    // atom is the sum of the shares of the blocks which carry it.

    std::vector<double> shares(natoms, 0.0);

    for (const auto &block : pattern.blocks())
    {
        const auto &c_atoms = block.c_atoms();

        if (c_atoms.empty()) continue;

        const auto share = static_cast<double>(block.number_of_elements() * sizeof(double)) /
                           static_cast<double>(c_atoms.size());

        for (const auto atom : c_atoms) shares[static_cast<size_t>(atom)] += share;
    }

    // NOTE: the atoms are gathered in the order they are given until the integrals
    // of a part reach the budget. An atom whose own integrals are above it is a part
    // of its own, as there is nothing smaller to divide.

    const CSimdThreeCenterElectronRepulsionDriver eri_drv;

    std::vector<CTripleSparsityPattern> parts;

    std::vector<int> atoms;

    double memory = 0.0;

    for (size_t atom = 0; atom < natoms; atom++)
    {
        const auto share = shares[atom];

        if (share <= 0.0) continue;

        if ((!atoms.empty()) && ((memory + share) > static_cast<double>(_budget / 2)))
        {
            parts.push_back(eri_drv.make_pattern(molecule, basis, aux_basis, threshold, atoms));

            atoms.clear();

            memory = 0.0;
        }

        atoms.push_back(static_cast<int>(atom));

        memory += share;
    }

    if (!atoms.empty()) parts.push_back(eri_drv.make_pattern(molecule, basis, aux_basis, threshold, atoms));

    return parts;
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
