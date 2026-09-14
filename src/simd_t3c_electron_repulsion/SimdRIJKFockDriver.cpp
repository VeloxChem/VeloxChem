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
#include <numeric>
#include <stdexcept>
#include <string>

#include "DenseIndexFunc.hpp"
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

/// @brief Writes the phases of a build of the direct mode, and their share of it.
/// @param times The times gathered over the calls the build was made of.
static auto
report_direct(const CDirectTimes &times) -> void
{
    static size_t builds = 0;

    builds++;

    const auto accounted = times.allocate + times.integrals_a + times.transform + times.closure +
                           times.copies + times.solve + times.exchange + times.integrals_b + times.coulomb;

    const char *names[] = {"allocate", "integrals a", "transform",   "closure", "copies",
                           "solve",    "exchange",    "integrals b", "coulomb", "rest"};

    const double values[] = {times.allocate, times.integrals_a, times.transform,   times.closure,
                             times.copies,   times.solve,       times.exchange,    times.integrals_b,
                             times.coulomb,  times.total - accounted};

    std::printf("RIJK build %zu on %d threads, %.3f s\n", builds, omp::get_number_of_threads(), times.total);

    for (size_t i = 0; i < 10; i++)
    {
        std::printf("RIJK   %-12s %9.3f s %6.1f %%\n", names[i], values[i],
                    (times.total > 0.0) ? 100.0 * values[i] / times.total : 0.0);
    }

    std::fflush(stdout);
}

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
CSimdRIJKFockDriver::required_memory(const CMolecule        &molecule,
                                     const CMolecularBasis  &basis,
                                     const CMolecularBasis  &aux_basis,
                                     const double            threshold,
                                     const std::vector<int> &aux_atoms) const -> size_t
{
    // NOTE: the pattern of the B vectors is the pattern of the three-center
    // integrals, as the metric is dense and the transformation of the auxiliary
    // side keeps every atom which survives.

    // NOTE: the memory answered is the memory of the atoms asked for, which is the
    // memory of this rank when the auxiliary basis is divided over a communicator.
    // Answering the memory of the whole molecule there would put every rank on the
    // direct way for a calculation each of them holds a fitting share of.

    const CSimdThreeCenterElectronRepulsionDriver eri_drv;

    const auto pattern = aux_atoms.empty() ? eri_drv.make_pattern(molecule, basis, aux_basis, threshold)
                                           : eri_drv.make_pattern(molecule, basis, aux_basis, threshold, aux_atoms);

    size_t nvalues = 0;

    for (size_t i = 0; i < static_cast<size_t>(pattern.number_of_blocks()); i++)
    {
        nvalues += pattern.block(i).number_of_elements();
    }

    return nvalues * sizeof(double);
}

namespace {

/// @brief Measures what each atom of the auxiliary side of a pattern carries.
/// @param pattern The sparsity pattern to measure.
/// @param natoms The number of atoms of the molecule.
/// @return The memory of the values of each atom, in bytes.
/// @note A block holds as many values for one of its atoms on the auxiliary side as
/// for any other, so its memory divides evenly over them and the memory of an atom is
/// the sum of the shares of the blocks which carry it.
static auto
atom_shares(const CTripleSparsityPattern &pattern, const size_t natoms) -> std::vector<double>
{
    std::vector<double> shares(natoms, 0.0);

    for (const auto &block : pattern.blocks())
    {
        const auto &c_atoms = block.c_atoms();

        if (c_atoms.empty()) continue;

        const auto share =
            static_cast<double>(block.number_of_elements() * sizeof(double)) / static_cast<double>(c_atoms.size());

        for (const auto atom : c_atoms) shares[static_cast<size_t>(atom)] += share;
    }

    return shares;
}

/// @brief Gets the dense indices of the auxiliary basis functions of given atoms.
/// @param aux_basis The auxiliary molecular basis.
/// @param atoms The atoms, as their indices in the molecule, or none of them for all
/// of them.
/// @return The indices, in ascending order.
/// @note The dense index runs over the angular momenta of the whole molecule before
/// it runs over the atoms, so the functions of one atom are scattered through it and
/// the functions of a set of atoms are a set rather than a range.
static auto
aux_functions_of(const CMolecularBasis &aux_basis, const std::vector<int> &atoms) -> std::vector<size_t>
{
    const auto set_indices = aux_basis.basis_sets_indices();

    const auto natoms = set_indices.size();

    const auto indices = denseidx::index_functions(aux_basis);

    const auto starts = denseidx::make_dense_starts(aux_basis);

    const auto strides = denseidx::make_dense_strides(aux_basis);

    const auto nmoms = static_cast<size_t>(aux_basis.max_angular_momentum() + 1);

    std::vector<int> all_atoms;

    if (atoms.empty())
    {
        all_atoms.reserve(natoms);

        for (size_t atom = 0; atom < natoms; atom++) all_atoms.push_back(static_cast<int>(atom));
    }

    std::vector<size_t> functions;

    for (const auto atom : (atoms.empty() ? all_atoms : atoms))
    {
        const auto index = static_cast<size_t>(atom);

        for (const auto [lc, kc] : indices[static_cast<size_t>(set_indices[index])])
        {
            const auto lval = static_cast<size_t>(lc);

            for (size_t mc = 0; mc < static_cast<size_t>(2 * lc + 1); mc++)
            {
                functions.push_back(starts[index * nmoms + lval] + kc + mc * strides[lval]);
            }
        }
    }

    std::sort(functions.begin(), functions.end());

    return functions;
}

/// @brief Forms the metric a way of building asks for, and the way it is for.
/// @param molecule The molecule to compute the metric of.
/// @param aux_basis The auxiliary molecular basis.
/// @param metric_threshold The threshold below which a direction is dropped.
/// @param use_inverse_square_root Whether to invert the square root of the metric.
/// @param mode The way of building the metric is for.
/// @param two_center_time Where to add the time of the two-center integrals, if
/// anywhere.
/// @param metric_time Where to add the time of the inversion, if anywhere.
/// @return The metric, and the way it is for.
static auto
form_metric(const CMolecule       &molecule,
            const CMolecularBasis &aux_basis,
            const double           metric_threshold,
            const bool             use_inverse_square_root,
            const rimode           mode,
            double                *two_center_time,
            double                *metric_time) -> std::pair<CPackedMatrix, rimode>
{
    const auto mark_two_center = prof_clock::now();

    const auto two_center = CSimdTwoCenterElectronRepulsionDriver().compute(molecule, aux_basis);

    if (two_center_time) *two_center_time += prof_since(mark_two_center);

    const auto mark_metric = prof_clock::now();

    // NOTE: both forms of the metric close the resolution of the identity, and the
    // Cholesky factor costs an order of magnitude less, so it is tried first. A
    // fitting basis which is close to linearly dependent has none, and the square
    // root is inverted in its place, dropping the directions which carry nothing.

    // NOTE: the direct way solves with the factor rather than multiplying by its
    // inverse, so it is the factor which is kept. The inverted square root is
    // taken for a metric which has no factor, in either way, as the B vectors
    // formed with it close the same sum.

    if (mode == rimode::direct)
    {
        if (!use_inverse_square_root)
        {
            try
            {
                auto factor = packlin::cholesky_factor(two_center);

                if (metric_time) *metric_time += prof_since(mark_metric);

                return {std::move(factor), rimode::direct};
            }
            catch (const std::runtime_error &)
            {
                errors::msg(std::string("RIJKFockDriver: The metric of the fitting basis has no Cholesky factor. The "
                                        "direct way needs one, so the B vectors are formed with the inverted square "
                                        "root of the metric instead and the calculation is held in memory."),
                            "Warning");
            }
        }
        else
        {
            errors::msg(std::string("RIJKFockDriver: The direct way solves with the Cholesky factor of the metric, "
                                    "which the inverted square root is not, so the calculation is held in memory."),
                        "Warning");
        }
    }

    auto metric = CPackedMatrix();

    if (use_inverse_square_root)
    {
        metric = packlin::inverse_square_root(two_center, metric_threshold);
    }
    else
    {
        try
        {
            metric = packlin::cholesky_inverse(two_center);
        }
        catch (const std::runtime_error &)
        {
            errors::msg(std::string("RIJKFockDriver: The metric of the fitting basis has no Cholesky factor, so its "
                                    "square root is inverted instead. This is a nearly linearly dependent fitting "
                                    "basis."),
                        "Warning");

            metric = packlin::inverse_square_root(two_center, metric_threshold);
        }
    }

    if (metric_time) *metric_time += prof_since(mark_metric);

    return {std::move(metric), rimode::in_memory};
}

}  // namespace

auto
CSimdRIJKFockDriver::make_metric(const CMolecule       &molecule,
                                 const CMolecularBasis &aux_basis,
                                 const double           metric_threshold,
                                 const bool             use_inverse_square_root,
                                 const rimode           mode) const -> std::pair<CPackedMatrix, rimode>
{
    errors::assertMsgCritical(mode != rimode::automatic,
                              std::string("RIJKFockDriver: The metric is formed for a named way of building"));

    return form_metric(molecule, aux_basis, metric_threshold, use_inverse_square_root, mode, nullptr, nullptr);
}

auto
CSimdRIJKFockDriver::prepare(const CMolecule       &molecule,
                             const CMolecularBasis &basis,
                             const CMolecularBasis &aux_basis,
                             const double           threshold,
                             const size_t           memory_budget,
                             const double           metric_threshold,
                             const bool             use_inverse_square_root,
                             const rimode           mode,
                             const std::vector<int> &aux_atoms,
                             const CPackedMatrix   &metric,
                             const size_t           min_parts) -> void
{
    CPrepareProfile profile;

    const auto profile_start = prof_clock::now();

    const auto memory = required_memory(molecule, basis, aux_basis, threshold, aux_atoms);

    _budget = memory_budget;

    // NOTE: the memory is answered from the sparsity pattern, before any integral
    // is computed, so a molecule whose B vectors do not fit is put on the direct
    // way at once rather than after the work of forming them.

    _mode = (mode == rimode::automatic) ? ((memory > memory_budget) ? rimode::direct : rimode::in_memory) : mode;

    _molecule = molecule;

    _basis = basis;

    _aux_basis = aux_basis;

    // NOTE: a metric given by the caller is taken as it is, and the way of building
    // with it must be named, as the fallbacks which change that way have already
    // been taken where it was formed. This is what lets the ranks of a communicator
    // share one metric: the master forms it with make_metric, which answers the way
    // as well, and hands both to every rank rather than each of them inverting the
    // same matrix and racing for the same fallback.

    const auto given = (metric.number_of_elements() > 0);

    errors::assertMsgCritical(!(given && (_mode == rimode::automatic)),
                              std::string("RIJKFockDriver: A metric given to the driver is for a named way of building"));

    auto [formed, formed_mode] =
        given ? std::pair<CPackedMatrix, rimode>{metric, _mode}
              : form_metric(molecule, aux_basis, metric_threshold, use_inverse_square_root, _mode, &profile.two_center, &profile.metric);

    _mode = formed_mode;

    if (_mode == rimode::direct)
    {
        // NOTE: the direct way divides the auxiliary basis into parts to hold the
        // half transformed integrals, but every part adds into the same array and
        // one triangular solve over the whole auxiliary basis follows, so a part is
        // not a summand: a rank given some of the atoms would solve a forward
        // substitution missing the rows above its own and answer a Fock matrix
        // which is not a share of anything. Dividing this way over a communicator
        // is refused here rather than silently answered wrongly.

        errors::assertMsgCritical(aux_atoms.empty(),
                                  std::string("RIJKFockDriver: The direct way cannot be divided over the atoms of the "
                                              "auxiliary basis, as its triangular solve reaches across all of them"));

        const auto mark_pattern = prof_clock::now();

        const CSimdThreeCenterElectronRepulsionDriver eri_drv;

        const auto pattern = eri_drv.make_pattern(molecule, basis, aux_basis, threshold);

        _parts = _make_parts(molecule, basis, aux_basis, threshold, pattern, min_parts);

        profile.pattern += prof_since(mark_pattern);

        _factor = std::move(formed);

        // NOTE: the direct way holds no B vectors and refuses a division of the
        // atoms, so it sweeps the whole auxiliary basis and says so.

        _aux_functions = aux_functions_of(aux_basis, {});

        _bq_vectors = CSparseTensor();

        _w_vectors.clear();

        _prepared = true;

        profile.total = prof_since(profile_start);

        if (prof_wanted()) profile.report("direct");

        return;
    }

    _metric = std::move(formed);

    _parts.clear();

    _aux_functions = aux_functions_of(aux_basis, aux_atoms);

    const auto mark_bq_vectors = prof_clock::now();

    _bq_vectors = _drv.compute_bq_vectors(molecule, basis, aux_basis, _metric, threshold, aux_atoms);

    profile.bq_vectors += prof_since(mark_bq_vectors);

    _w_vectors.clear();

    _drv.release_squares();

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

    if (_mode == rimode::direct) return _compute_direct(coefficients, exchange_scaling_factor);

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

    // NOTE: the bound is the smaller of the constant and a share of the budget of
    // this driver, which is the share of a rank of whatever machine it is on. The
    // constant alone is a bound on a process, and a node given to several of them
    // would be promised its memory several times over.

    const auto allowance = std::min(_w_batch_memory, std::max(_budget / _w_batch_divisor, size_t{1}));

    const auto by_memory = std::max(size_t{1}, allowance / std::max(per_function, size_t{1}));

    const auto nheld = _aux_functions.size();

    const auto nbatch = std::min(nheld, std::max(_w_batch, by_memory));

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

    for (size_t first = 0; first < nheld; first += nbatch)
    {
        const auto last = std::min(first + nbatch, nheld);

        const auto count = last - first;

        // NOTE: the functions of this range of the set this driver holds, which is
        // what the transformation is asked for. A rank which holds a share of the
        // atoms therefore forms a share of the W matrices, and its exchange is the
        // square of that share rather than of the whole auxiliary basis with the
        // other ranks' functions in it as zeros.

        const auto batch = std::vector<size_t>(_aux_functions.begin() + static_cast<long>(first),
                                               _aux_functions.begin() + static_cast<long>(last));

        // NOTE: the last range is shorter than the others, and the storage is
        // handed to the transformation as the range it is asked to fill.

        const auto mark_transform = prof_clock::now();

        if (count == nbatch)
        {
            _drv.compute_w_vectors(_bq_vectors, _basis, _aux_basis, coefficients, batch, _w_vectors);

            profile.transform += prof_since(mark_transform);

            const auto mark_exchange = prof_clock::now();

            _drv.compute_exchange_matrix(_w_vectors, fock, -exchange_scaling_factor);

            profile.exchange += prof_since(mark_exchange);
        }
        else
        {
            auto tail = std::vector<CPackedMatrix>(_w_vectors.begin(), _w_vectors.begin() + static_cast<long>(count));

            _drv.compute_w_vectors(_bq_vectors, _basis, _aux_basis, coefficients, batch, tail);

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
CSimdRIJKFockDriver::compute_exchange(const CPackedMatrix &coefficients,
                                      const double         exchange_scaling_factor,
                                      const size_t         ofirst,
                                      const size_t         olast) -> std::pair<CPackedMatrix, std::vector<double>>
{
    errors::assertMsgCritical(_prepared, std::string("RIJKFockDriver: The driver has not been prepared"));

    errors::assertMsgCritical(_mode == rimode::direct,
                              std::string("RIJKFockDriver: The exchange pass belongs to the direct way, and the way "
                                          "which holds the B vectors forms its Fock matrices in one call"));

    const auto nao = _basis.dimensions_of_basis();

    const auto naux = _aux_basis.dimensions_of_basis();

    const auto norbitals = coefficients.number_of_columns();

    errors::assertMsgCritical((coefficients.get_type() == mat_t::general) && (coefficients.number_of_rows() == nao),
                              std::string("RIJKFockDriver: The orbital coefficients do not match the molecular basis"));

    errors::assertMsgCritical((ofirst <= olast) && (olast <= norbitals),
                              std::string("RIJKFockDriver: The range of orbitals is not a range of the coefficients"));

    // NOTE: the times of the whole build are gathered from here, as this is the
    // first of the calls a build is made of.

    _direct_times = CDirectTimes();

    const auto profile_start = prof_clock::now();

    auto fock = CPackedMatrix(nao, nao, mat_t::symmetric);

    fock.zero();

    std::vector<double> gamma(naux, 0.0);

    if (ofirst == olast)
    {
        _direct_times.total += prof_since(profile_start);

        return {std::move(fock), std::move(gamma)};
    }

    // NOTE: the half transformed integrals are the auxiliary basis by the basis
    // functions by the orbitals of a batch, and are held twice: once as the
    // matrices the transformation writes and once as the square the triangular
    // solve reads. The batch of orbitals is chosen so that the two fit.

    const auto per_orbital = 2 * naux * nao * sizeof(double);

    const auto ntaken = olast - ofirst;

    const auto nbatch = std::max(size_t{1}, std::min(ntaken, (_budget / 2) / std::max(per_orbital, size_t{1})));

    const auto *cvalues = coefficients.data();

    CSimdThreeCenterElectronRepulsionDriver eri_drv;

    // one sweep of the integrals for every batch of the orbitals asked for

    for (size_t first = ofirst; first < olast; first += nbatch)
    {
        const auto last = std::min(first + nbatch, olast);

        const auto ncols = last - first;

        auto batch = CPackedMatrix(nao, ncols, mat_t::general);

        auto *bvalues = batch.data();

        for (size_t irow = 0; irow < nao; irow++)
        {
            std::copy(cvalues + irow * norbitals + first, cvalues + irow * norbitals + last,
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

        _direct_times.allocate += prof_since(mark_allocate);

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

            _direct_times.integrals_a += prof_since(mark_integrals);

            const auto mark_transform = prof_clock::now();

            _drv.compute_w_vectors(integrals, _basis, _aux_basis, batch, 0, naux, half, true);

            _direct_times.transform += prof_since(mark_transform);
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

        _direct_times.closure += prof_since(mark_closure);

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

        _direct_times.copies += prof_since(mark_copies);

        const auto mark_solve = prof_clock::now();

        _solve_factor(stacked.get(), naux, width, false);

        _direct_times.solve += prof_since(mark_solve);

        mark_copies = prof_clock::now();

#pragma omp parallel for schedule(static) if (ncopies > 1)
        for (int iq = 0; iq < ncopies; iq++)
        {
            const auto q = static_cast<size_t>(iq);

            std::copy(stacked.get() + q * width, stacked.get() + (q + 1) * width, half[q].data());
        }

        _direct_times.copies += prof_since(mark_copies);

        const auto mark_exchange = prof_clock::now();

        _drv.compute_exchange_matrix(half, fock, -exchange_scaling_factor);

        _direct_times.exchange += prof_since(mark_exchange);
    }

    _direct_times.total += prof_since(profile_start);

    return {std::move(fock), std::move(gamma)};
}

auto
CSimdRIJKFockDriver::solve_fitting(std::vector<double> gamma) -> std::vector<double>
{
    errors::assertMsgCritical(_prepared, std::string("RIJKFockDriver: The driver has not been prepared"));

    errors::assertMsgCritical(_mode == rimode::direct,
                              std::string("RIJKFockDriver: The fitting belongs to the direct way"));

    const auto naux = _aux_basis.dimensions_of_basis();

    errors::assertMsgCritical(gamma.size() == naux,
                              std::string("RIJKFockDriver: The right hand side of the fitting is not one value per "
                                          "auxiliary basis function"));

    const auto profile_start = prof_clock::now();

    // the coefficients of the fitting, from the factor and its transpose

    _solve_factor(gamma.data(), naux, 1, false);

    _solve_factor(gamma.data(), naux, 1, true);

    _direct_times.solve += prof_since(profile_start);

    _direct_times.total += prof_since(profile_start);

    return gamma;
}

auto
CSimdRIJKFockDriver::compute_coulomb(const std::vector<double> &gamma,
                                     const std::vector<int>    &parts,
                                     CPackedMatrix             &matrix) -> void
{
    errors::assertMsgCritical(_prepared, std::string("RIJKFockDriver: The driver has not been prepared"));

    errors::assertMsgCritical(_mode == rimode::direct,
                              std::string("RIJKFockDriver: The Coulomb pass belongs to the direct way"));

    const auto nao = _basis.dimensions_of_basis();

    errors::assertMsgCritical(gamma.size() == _aux_basis.dimensions_of_basis(),
                              std::string("RIJKFockDriver: The coefficients of the fitting are not one value per "
                                          "auxiliary basis function"));

    errors::assertMsgCritical((matrix.number_of_rows() == nao) && (matrix.number_of_columns() == nao),
                              std::string("RIJKFockDriver: The matrix to add the Coulomb matrix to does not match the "
                                          "molecular basis"));

    const auto profile_start = prof_clock::now();

    CSimdThreeCenterElectronRepulsionDriver eri_drv;

    for (const auto index : parts)
    {
        errors::assertMsgCritical((index >= 0) && (static_cast<size_t>(index) < _parts.size()),
                                  std::string("RIJKFockDriver: The Coulomb matrix was asked for a part of the "
                                              "auxiliary basis which the driver does not sweep"));

        const auto &pattern = _parts[static_cast<size_t>(index)];

        auto integrals = CSparseTensor(pattern);

        integrals.allocate();

        auto distributor = CSimdT3CDistributor<CSparseTensor>(&integrals);

        const auto mark_integrals = prof_clock::now();

        eri_drv.compute(pattern, _molecule, _basis, _aux_basis, distributor);

        _direct_times.integrals_b += prof_since(mark_integrals);

        const auto mark_coulomb = prof_clock::now();

        // NOTE: the blocks of atom pairs of one batch write elements of the matrix
        // no other batch writes, so the parts are added and nothing is counted
        // twice. The Coulomb matrix enters twice, as the density is of one spin.

        auto part_fock = _drv.compute_fock_matrix(integrals, _basis, _aux_basis, gamma);

        auto       *values = matrix.data();

        const auto *added = part_fock.data();

        const auto nvalues = static_cast<int>(matrix.number_of_elements());

#pragma omp parallel for schedule(static) if (nvalues > 1)
        for (int i = 0; i < nvalues; i++)
        {
            values[static_cast<size_t>(i)] += 2.0 * added[static_cast<size_t>(i)];
        }

        _direct_times.coulomb += prof_since(mark_coulomb);
    }

    // NOTE: the times of the whole build are reported from here, as this is the
    // last of the calls a build is made of.

    _direct_times.total += prof_since(profile_start);

    if (prof_wanted()) report_direct(_direct_times);
}

auto
CSimdRIJKFockDriver::aux_atom_weights(const CMolecule       &molecule,
                                      const CMolecularBasis &basis,
                                      const CMolecularBasis &aux_basis,
                                      const double           threshold) const -> std::vector<double>
{
    const auto pattern = CSimdThreeCenterElectronRepulsionDriver().make_pattern(molecule, basis, aux_basis, threshold);

    return atom_shares(pattern, static_cast<size_t>(molecule.number_of_atoms()));
}

auto
CSimdRIJKFockDriver::number_of_aux_functions() const -> size_t
{
    return _aux_functions.size();
}

auto
CSimdRIJKFockDriver::number_of_parts() const -> size_t
{
    return _parts.size();
}

auto
CSimdRIJKFockDriver::_compute_direct(const CPackedMatrix &coefficients,
                                     const double         exchange_scaling_factor) -> CPackedMatrix
{
    // NOTE: the three calls of a build, in a row, over the whole of the orbitals
    // and the whole of the auxiliary basis. A caller dividing the work over a
    // communicator makes the same three calls with a range and a share of the
    // parts each, and adds the matrices; this is here so that one rank is one
    // call, and so that the two ways read the same from the outside.

    auto [fock, gamma] = compute_exchange(coefficients, exchange_scaling_factor, 0,
                                          coefficients.number_of_columns());

    gamma = solve_fitting(std::move(gamma));

    std::vector<int> parts(_parts.size());

    std::iota(parts.begin(), parts.end(), 0);

    compute_coulomb(gamma, parts, fock);

    return fock;
}

auto
CSimdRIJKFockDriver::_make_parts(const CMolecule              &molecule,
                                 const CMolecularBasis        &basis,
                                 const CMolecularBasis        &aux_basis,
                                 const double                  threshold,
                                 const CTripleSparsityPattern &pattern,
                                 const size_t                  min_parts) const -> std::vector<CTripleSparsityPattern>
{
    const auto natoms = static_cast<size_t>(molecule.number_of_atoms());

    const auto shares = atom_shares(pattern, natoms);

    // NOTE: the parts are cut at whichever is the smaller of what the memory allows
    // and an equal division into the number asked for. A machine with memory to spare
    // gives one part, and the caller which divides the Coulomb pass over the ranks of
    // a communicator then has one part for all of them: one rank sweeps the integrals
    // a second time and the others wait. Cutting finer costs nothing, as the parts are
    // a division of the same atoms and their integrals are the same integrals however
    // they are grouped.

    const auto total = std::accumulate(shares.begin(), shares.end(), 0.0);

    const auto by_memory = static_cast<double>(_budget / 2);

    const auto by_parts = total / static_cast<double>(std::max(min_parts, size_t{1}));

    const auto cut = std::min(by_memory, by_parts);

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

        if ((!atoms.empty()) && ((memory + share) > cut))
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
