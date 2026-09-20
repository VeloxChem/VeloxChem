//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//

#include "SimdRIFockCommon.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>

#include "DenseIndexFunc.hpp"
#include "ErrorHandler.hpp"
#include "Eigen/Dense"
#include "MathLibrary.hpp"
#include "PackedLinearAlgebra.hpp"
#include "SimdTwoCenterElectronRepulsionDriver.hpp"
#include "StringFormat.hpp"

namespace simdri {

namespace {

using prof_clock = std::chrono::steady_clock;

/// @brief The seconds since a mark.
auto
prof_since(const prof_clock::time_point &mark) -> double
{
    return std::chrono::duration<double>(prof_clock::now() - mark).count();
}

}  // namespace

/// @brief Measures what each atom of the auxiliary side of a pattern carries.
/// @param pattern The sparsity pattern to measure.
/// @param natoms The number of atoms of the molecule.
/// @return The memory of the values of each atom, in bytes.
/// @note A block holds as many values for one of its atoms on the auxiliary side as
/// for any other, so its memory divides evenly over them and the memory of an atom is
/// the sum of the shares of the blocks which carry it.
auto
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
auto
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

/// @brief Inverts a metric by the route asked for, falling back where there is no
/// Cholesky factor to be had.
/// @param two_center The metric of the fitting basis.
/// @param metric_threshold The threshold below which a direction is dropped.
/// @param use_inverse_square_root Whether to invert the square root of the metric.
/// @param what What the metric is, named in the warning a fallback prints, so that a
/// caller inverting two of them says which of the two it was.
/// @return The inverted metric.
auto
invert_metric(const CPackedMatrix &two_center,
              const double         metric_threshold,
              const bool           use_inverse_square_root,
              const std::string   &what) -> CPackedMatrix
{
    if (use_inverse_square_root) return packlin::inverse_square_root(two_center, metric_threshold);

    try
    {
        return packlin::cholesky_inverse(two_center);
    }
    catch (const std::runtime_error &)
    {
        errors::msg(std::string("RIJKFockDriver: The ") + what +
                        std::string(" has no Cholesky factor, so its square root is inverted instead. This is a "
                                    "nearly linearly dependent fitting basis."),
                    "Warning");

        return packlin::inverse_square_root(two_center, metric_threshold);
    }
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
auto
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
                errors::msg(std::string("RIJKFockDriver: The metric of the fitting basis has no Cholesky factor, so "
                                        "its square root is inverted instead and multiplied by. This is a nearly "
                                        "linearly dependent fitting basis."),
                            "Warning");
            }
        }

        // NOTE: the direct way solves the Cholesky factor against the half
        // transformed integrals where it has one, and multiplies by the inverted
        // square root where that is what was asked for, or where there is no factor
        // to be had. The two close the same sum: solving the factor gives B with
        // B^T B equal to A^T V^-1 A, and so does multiplying by the root, the root
        // being its own transpose. The root costs twice the arithmetic and is a
        // product rather than a substitution, which is the trade. Falling back to
        // holding the B vectors instead, which is what this did, asks for the
        // memory the direct way was chosen for want of.

        auto root = packlin::inverse_square_root(two_center, metric_threshold);

        if (metric_time) *metric_time += prof_since(mark_metric);

        return {std::move(root), rimode::direct};
    }

    auto metric = invert_metric(two_center, metric_threshold, use_inverse_square_root,
                                "metric of the fitting basis");

    if (metric_time) *metric_time += prof_since(mark_metric);

    return {std::move(metric), rimode::in_memory};
}

}  // namespace simdri
