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
#include "ScreeningFunc.hpp"
#include "SimdThreeCenterElectronRepulsionDriver.hpp"
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
                             const bool             use_inverse_square_root) -> void
{
    const auto memory = required_memory(molecule, basis, aux_basis, threshold);

    // NOTE: the check is made before the integrals are computed, so that a
    // calculation which does not fit is told the numbers rather than reaching the
    // allocator after the work of forming them.

    if (memory > memory_budget)
    {
        const auto needed = static_cast<double>(memory) / (1024.0 * 1024.0 * 1024.0);

        const auto budget = static_cast<double>(memory_budget) / (1024.0 * 1024.0 * 1024.0);

        // NOTE: this is an exception rather than a critical error. The check exists
        // so that a caller can choose another basis or another machine, and a
        // critical error would end the interpreter instead of letting it choose.

        throw std::runtime_error(std::string("RIJKFockDriver.prepare: The B vectors need ") + std::to_string(needed) +
                                 std::string(" GB and the budget is ") + std::to_string(budget) + std::string(" GB"));
    }

    _basis = basis;

    _aux_basis = aux_basis;

    const auto two_center = CSimdTwoCenterElectronRepulsionDriver().compute(molecule, aux_basis);

    // NOTE: both forms of the metric close the resolution of the identity, and the
    // Cholesky factor costs an order of magnitude less, so it is tried first. A
    // fitting basis which is close to linearly dependent has none, and the square
    // root is inverted in its place, dropping the directions which carry nothing.

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
CSimdRIJKFockDriver::get_bq_vectors() const -> const CSparseTensor &
{
    return _bq_vectors;
}

auto
CSimdRIJKFockDriver::get_metric() const -> const CPackedMatrix &
{
    return _metric;
}
