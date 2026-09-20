//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//

#include "SimdRIJFockDriver.hpp"

#include <algorithm>
#include <string>

#include "ErrorHandler.hpp"
#include "SimdThreeCenterElectronRepulsionDriver.hpp"
#include "SimdTwoCenterElectronRepulsionDriver.hpp"

auto
CSimdRIJFockDriver::required_memory(const CMolecule        &molecule,
                                    const CMolecularBasis  &basis,
                                    const CMolecularBasis  &aux_basis,
                                    const double            threshold,
                                    const std::vector<int> &aux_atoms) const -> size_t
{
    return simdri::pattern_memory(molecule, basis, aux_basis, threshold, aux_atoms);
}

auto
CSimdRIJFockDriver::make_metric(const CMolecule       &molecule,
                                const CMolecularBasis &aux_basis,
                                const double           metric_threshold) const -> CPackedMatrix
{
    const CSimdTwoCenterElectronRepulsionDriver t2c_drv;

    return simdri::invert_metric_full(t2c_drv.compute(molecule, aux_basis), metric_threshold);
}

auto
CSimdRIJFockDriver::prepare(const CMolecule       &molecule,
                            const CMolecularBasis &basis,
                            const CMolecularBasis &aux_basis,
                            const double           threshold,
                            const size_t           memory_budget,
                            const double           metric_threshold,
                            const rimode           mode,
                            const CPackedMatrix   &metric,
                            const size_t           nodes) -> void
{
    _molecule = molecule;

    _basis = basis;

    _aux_basis = aux_basis;

    _metric = (metric.number_of_rows() > 0) ? metric : make_metric(molecule, aux_basis, metric_threshold);

    const auto naux = aux_basis.dimensions_of_basis();

    errors::assertMsgCritical((_metric.number_of_rows() == naux) && (_metric.number_of_columns() == naux),
                              std::string("RIJFockDriver: The metric does not match the auxiliary basis"));

    // NOTE: the parts are what a sweep takes one at a time, and the ranks of a
    // communicator take some parts each, so there have to be at least as many parts
    // as there are ranks or a rank is given nothing to do.

    const CSimdThreeCenterElectronRepulsionDriver eri_drv;

    const auto pattern = eri_drv.make_pattern(molecule, basis, aux_basis, threshold);

    _parts = simdri::make_parts(molecule, basis, aux_basis, threshold, pattern, std::max(nodes, size_t{1}),
                                memory_budget);

    // NOTE: holding the integrals is a choice about memory and nothing else. Either
    // way the same sweeps run over the same parts; the way which holds them forms
    // them once and the way which does not forms them on every build. There is no
    // transformation between the two, which is why the choice is so much less
    // consequential here than it is for a driver which also forms the exchange.

    _mode = mode;

    if (_mode == rimode::automatic)
    {
        _mode = (simdri::pattern_memory(molecule, basis, aux_basis, threshold, {}) > memory_budget) ? rimode::direct
                                                                                                   : rimode::in_memory;
    }

    _integrals.clear();

    if (_mode == rimode::in_memory)
    {
        _integrals.reserve(_parts.size());

        for (const auto &part : _parts)
        {
            _integrals.push_back(simdri::integrals_of_part(part, _molecule, _basis, _aux_basis));
        }
    }

    _prepared = true;
}

auto
CSimdRIJFockDriver::is_prepared() const -> bool
{
    return _prepared;
}

auto
CSimdRIJFockDriver::get_mode() const -> rimode
{
    return _mode;
}

auto
CSimdRIJFockDriver::get_metric() const -> const CPackedMatrix &
{
    return _metric;
}

auto
CSimdRIJFockDriver::number_of_parts() const -> size_t
{
    return _parts.size();
}

auto
CSimdRIJFockDriver::_check_part(const int index) const -> void
{
    errors::assertMsgCritical((index >= 0) && (static_cast<size_t>(index) < _parts.size()),
                              std::string("RIJFockDriver: A part of the auxiliary basis was asked for which the "
                                          "driver does not sweep"));
}

auto
CSimdRIJFockDriver::_part_integrals(const size_t index) -> CSparseTensor
{
    // NOTE: the way which holds them hands back a copy rather than a reference. The
    // sweeps do not write to the integrals, so a reference would do and would save
    // the copy; it is a copy here so that the two ways have one signature between
    // them and the sweeps below read the same.

    if (_mode == rimode::in_memory) return _integrals[index];

    return simdri::integrals_of_part(_parts[index], _molecule, _basis, _aux_basis);
}

auto
CSimdRIJFockDriver::compute_gamma(const CPackedMatrix &density, const std::vector<int> &parts) -> std::vector<double>
{
    errors::assertMsgCritical(_prepared, std::string("RIJFockDriver: The driver has not been prepared"));

    const auto nao = _basis.dimensions_of_basis();

    errors::assertMsgCritical((density.number_of_rows() == nao) && (density.number_of_columns() == nao),
                              std::string("RIJFockDriver: The density does not match the molecular basis"));

    auto gamma = std::vector<double>(_aux_basis.dimensions_of_basis(), 0.0);

    for (const auto index : parts)
    {
        _check_part(index);

        const auto integrals = _part_integrals(static_cast<size_t>(index));

        // NOTE: a part carries the auxiliary functions of its own atoms and nothing
        // else, so what comes back is zero everywhere but there and the parts add
        // into one vector without any of them overwriting another.

        const auto partial = _drv.compute_y_vector(integrals, _basis, _aux_basis, density);

        errors::assertMsgCritical(partial.size() == gamma.size(),
                                  std::string("RIJFockDriver: A part answered a right hand side which is not one "
                                              "value per auxiliary basis function"));

        for (size_t i = 0; i < gamma.size(); i++) gamma[i] += partial[i];
    }

    return gamma;
}

auto
CSimdRIJFockDriver::solve_fitting(const std::vector<double> &gamma) const -> std::vector<double>
{
    errors::assertMsgCritical(_prepared, std::string("RIJFockDriver: The driver has not been prepared"));

    const auto naux = _aux_basis.dimensions_of_basis();

    errors::assertMsgCritical(gamma.size() == naux,
                              std::string("RIJFockDriver: The right hand side of the fitting is not one value per "
                                          "auxiliary basis function"));

    // NOTE: the metric is held in the packed format and the multiply wants a square,
    // so it is expanded here. It is the square of the auxiliary basis, which is small
    // beside the integrals of any part of it.

    auto dense = std::vector<double>(naux * naux, 0.0);

    _metric.to_dense(dense.data());

    auto fitted = std::vector<double>(naux, 0.0);

    for (size_t i = 0; i < naux; i++)
    {
        double sum = 0.0;

        const auto *row = dense.data() + i * naux;

        for (size_t j = 0; j < naux; j++) sum += row[j] * gamma[j];

        fitted[i] = sum;
    }

    return fitted;
}

auto
CSimdRIJFockDriver::compute_coulomb(const std::vector<double> &gamma,
                                    const std::vector<int>    &parts,
                                    CPackedMatrix             &matrix) -> void
{
    errors::assertMsgCritical(_prepared, std::string("RIJFockDriver: The driver has not been prepared"));

    const auto nao = _basis.dimensions_of_basis();

    errors::assertMsgCritical(gamma.size() == _aux_basis.dimensions_of_basis(),
                              std::string("RIJFockDriver: The coefficients of the fitting are not one value per "
                                          "auxiliary basis function"));

    errors::assertMsgCritical((matrix.number_of_rows() == nao) && (matrix.number_of_columns() == nao),
                              std::string("RIJFockDriver: The matrix to add the Coulomb matrix to does not match the "
                                          "molecular basis"));

    for (const auto index : parts)
    {
        _check_part(index);

        const auto integrals = _part_integrals(static_cast<size_t>(index));

        // NOTE: the blocks of atom pairs of one part write elements no other part
        // writes, so the parts are added and nothing is counted twice.

        const auto part_fock = _drv.compute_fock_matrix(integrals, _basis, _aux_basis, gamma);

        auto *values = matrix.data();

        const auto *added = part_fock.data();

        const auto nvalues = static_cast<int>(matrix.number_of_elements());

#pragma omp parallel for schedule(static) if (nvalues > 1)
        for (int i = 0; i < nvalues; i++)
        {
            values[static_cast<size_t>(i)] += added[static_cast<size_t>(i)];
        }
    }
}

auto
CSimdRIJFockDriver::compute(const CPackedMatrix &density) -> CPackedMatrix
{
    errors::assertMsgCritical(_prepared, std::string("RIJFockDriver: The driver has not been prepared"));

    auto parts = std::vector<int>(_parts.size());

    for (size_t i = 0; i < parts.size(); i++) parts[i] = static_cast<int>(i);

    const auto gamma = solve_fitting(compute_gamma(density, parts));

    const auto nao = _basis.dimensions_of_basis();

    auto matrix = CPackedMatrix(nao, nao, mat_t::symmetric);

    matrix.zero();

    compute_coulomb(gamma, parts, matrix);

    return matrix;
}
